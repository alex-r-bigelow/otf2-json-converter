#!/usr/bin/env python3
import sys
import argparse
import subprocess
import re
import json
import heapq
import newick

parser = argparse.ArgumentParser(description='Collect stdout from a phylanx')
parser.add_argument('-i', '--input', dest='input', default=sys.stdin, type=argparse.FileType('r'), nargs='?',
                    help='stdout from phylanx as a file; alternatively, you can omit this argument and pipe the phylanx output directly into this script')
parser.add_argument('-o', '--otf2', dest='otf2', default='test_data/OTF2_archive/APEX.otf2',
                    help='The input otf2 trace file (e.g. test_data/OTF2_archive/APEX.otf2)')
parser.add_argument('-e', '--events', dest='events', action="store_true",
                    help='Include all discrete trace events (e.g. ENTER, LEAVE, MPI_SEND, etc. are separate)')

args = parser.parse_args()
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# Output the opening JSON character
sys.stdout.write('{')

# Parse the tree and performance data (waits for stdout to finish before attempting to parse the OTF2 trace)
treeModeParser = re.compile(r'Tree information for function:')
codeLinkedEventParser = re.compile(r'(/phylanx/[^\$]*)\$([^\$]*)\$([^\$]*)\$([^\$]*)\$([^\$]*)')

regions = {}
def dumpTree(node, spaces, parent=None):
    # Create the hashed region object
    regionName = node.name.strip()
    assert regionName not in regions
    regions[regionName] = {}
    if parent is not None:
        regions[regionName]['parent'] = parent
    codeLinkedEvent = codeLinkedEventParser.match(regionName)
    if codeLinkedEvent is not None:
        regions[regionName]['name'] = codeLinkedEvent[1]
        regions[regionName]['token1'] = codeLinkedEvent[2]
        regions[regionName]['token2'] = codeLinkedEvent[3]
        regions[regionName]['token3'] = codeLinkedEvent[4]
        regions[regionName]['token4'] = codeLinkedEvent[5]

    # Create the tree hierarchy
    sys.stdout.write('\n' + spaces + '"name": "' + regionName + '",')
    if len(node.descendants) > 0:
        sys.stdout.write('\n' + spaces + '"children": [')
        for index, child in enumerate(node.descendants):
            if index > 0:
                sys.stdout.write(',')
            sys.stdout.write('\n' + spaces + '  {')
            dumpTree(child, spaces + '    ', parent=regionName)
            sys.stdout.write('\n' + spaces + '  }')
        sys.stdout.write('\n' + spaces + ']')
    else:
        sys.stdout.write('\n' + spaces + '"children": []')

mode = None
for line in args.input:
    if mode is None:
        if treeModeParser.match(line):
            mode = 'tree'
    elif mode == 'tree':
        sys.stdout.write('\n  "tree": {')
        dumpTree(newick.loads(line)[0], '    ')
        sys.stdout.write('\n  }')
        mode = None
    else:
        assert False

# Parse the OTF2 trace, output non-ENTER/LEAVE events directly as we encounter them so they don't stick around in memory
otfPrint = subprocess.Popen(['otf2-print', args.otf2], stdout=subprocess.PIPE)

eventLineParser = re.compile(r'^(\S+)\s+(\d+)\s+(\d+)\s+(.*)$')
attrParsers = {
  'ENTER': r'(Region): "([^"]*)"',
  'LEAVE': r'(Region): "([^"]*)"',
  'METRIC': r'Value: \("([^"]*)" <\d+>; [^;]*; ([^\)]*)',
  'MPI_SEND': r'([^:]*): ([^,]*)',
  'MPI_RECV': r'([^:]*): ([^,]*)'
}
addAttrLineParser = re.compile(r'^\s+ADDITIONAL ATTRIBUTES: (.*)$')
addAttrSplitter = re.compile(r'\), \(')
addAttrParser = re.compile(r'\(?"([^"]*)" <\d+>; [^;]*; ([^\)]*)')

if args.events is True:
    sys.stdout.write(',\n  "events": {')
currentEvent = None
numEvents = 0
locations = {}

def dumpEvent(currentEvent, numEvents):
    if currentEvent is not None:
        numEvents += 1
        if numEvents % 10000 == 0:
            eprint('.', end=''),
        if numEvents % 100000 == 0:
            eprint('processed %i events' % numEvents)

        regionName = currentEvent['Region']
        if not regionName in regions:
            regions[regionName] = {}
            # TODO: figure out parent region based on GUIDs
            regions['name'] = regionName
        if not 'eventCount' in regions[regionName]:
            regions[regionName]['eventCount'] = 0
        regions[regionName]['eventCount'] += 1

        if args.events is not None:
            # Dump the event to the file
            if numEvents > 1:
                sys.stdout.write(',')
            sys.stdout.write('\n    ' + json.dumps(currentEvent))

        if currentEvent['Event'] == 'ENTER' or currentEvent['Event'] == 'LEAVE':
            # Add (and sort) the enter / leave event into its location list
            if not currentEvent['Location'] in locations:
                locations[currentEvent['Location']] = []
            heapq.heappush(locations[currentEvent['Location']], (currentEvent['Timestamp'], currentEvent))
    return numEvents

for line in otfPrint.stdout:
    line = line.decode()
    eventLineMatch = eventLineParser.match(line)
    addAttrLineMatch = addAttrLineParser.match(line)
    if currentEvent is None and eventLineMatch is None:
        # This is a blank / header line
        continue

    if eventLineMatch is not None:
        # This is the beginning of a new event; dump the previous one
        numEvents = dumpEvent(currentEvent, numEvents)
        currentEvent = {}
        currentEvent['Event'] = eventLineMatch.group(1)
        currentEvent['Location'] = int(eventLineMatch.group(2))
        currentEvent['Timestamp'] = int(eventLineMatch.group(3))
        attrs = eventLineMatch.group(4)
        assert currentEvent['Event'] in attrParsers
        for attrMatch in re.finditer(attrParsers[currentEvent['Event']], attrs):
            currentEvent[attrMatch.group(1)] = attrMatch.group(2)
    else:
        # This line contains additional event attributes
        assert currentEvent is not None and addAttrLineMatch is not None
        attrList = addAttrSplitter.split(addAttrLineMatch.group(1))
        for attrStr in attrList:
            attr = addAttrParser.match(attrStr)
            assert attr is not None
            currentEvent[attr.group(1)] = attr.group(2)
# The last event will never have had a chance to be dumped:
dumpEvent(currentEvent, numEvents)
if args.events is True:
    sys.stdout.write('\n  }')

# Combine the sorted enter / leave events into ranges
sys.stdout.write(',\n  "ranges":{')
numRanges = 0
for eventList in locations.values():
    lastEvent = None
    for timestamp, event in eventList:
        if event['Event'] == 'ENTER':
            # Start a range (don't output anything)
            assert lastEvent is None
            lastEvent = event
        elif event['Event'] == 'LEAVE':
            # Finish a range
            assert lastEvent is not None
            currentRange = {}
            for attr, value in event.items():
                # For now, we assume Event, Timestamp, and Region are the only things
                # that can change between an ENTER / LEAVE pair
                if attr != 'Event' and attr != 'Timestamp' and attr != 'Region':
                    assert event[attr] == lastEvent[attr]
                    currentRange[attr] = value
            currentRange['enter'] = {
                'Timestamp': lastEvent['Timestamp'],
                'Region': lastEvent['Region']
            }
            currentRange['leave'] = {
                'Timestamp': event['Timestamp'],
                'Region': event['Region']
            }

            # Output the finished range
            numRanges += 1
            if numRanges % 10000 == 0:
                eprint('.', end=''),
            if numRanges % 100000 == 0:
                eprint('processed %i ranges' % numRanges)
            if numRanges > 0:
                sys.stdout.write(',')
            sys.stdout.write('\n    ' + json.dumps(currentRange))
            lastEvent = None
    # Make sure there are no trailing ENTER events
    assert lastEvent is None
# Finish the ranges dict
sys.stdout.write('\n  }')


# Output the regions hash
sys.stdout.write(',\n  "regions": {')
for index, (key, value) in enumerate(regions.items()):
    if index > 0:
        sys.stdout.write(',')
    sys.stdout.write('\n    "' + key + '": ' + json.dumps(value))
sys.stdout.write('\n  }')

# Finish the JSON output
sys.stdout.write('\n}')