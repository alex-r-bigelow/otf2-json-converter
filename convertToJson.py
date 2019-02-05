#!/usr/bin/env python3
import sys
import argparse
import subprocess
import re
import json
import heapq
import newick

parser = argparse.ArgumentParser(description='Collect stdout and/or OTF2 trace data from a phylanx run')
parser.add_argument('-i', '--input', dest='input', default=sys.stdin, type=argparse.FileType('r'), nargs='?',
                    help='stdout from phylanx run as a file; alternatively, you can omit this argument and pipe the phylanx output directly into this script')
parser.add_argument('-o', '--otf2', dest='otf2',
                    help='The input otf2 trace file (e.g. test_data/OTF2_archive/APEX.otf2)')
parser.add_argument('-e', '--events', dest='events', action='store_true',
                    help='Include all discrete trace events (e.g. ENTER, LEAVE, MPI_SEND, etc. are separate)')
parser.add_argument('-t', '--tree', dest='tree', action='store_true',
                    help='Include the tree (separate from the implicit one in regions)')
parser.add_argument('-r', '--omit_ranges', dest='ranges', action='store_false',
                    help='Suppress trace ranges (ENTER and LEAVE are combined)')
parser.add_argument('-l', '--omit_links', dest='region_links', action='store_false',
                    help='Suppress the region links (separate from the implicit links in regions)')

args = parser.parse_args()

# Output the opening JSON character
sys.stdout.write('{')

# Helper tools for outputting / logging data
def log(value, end='\n'):
    sys.stderr.write(value + end)
    sys.stderr.flush()

firstRootElement = True
def writeRootComma(condition):
    global firstRootElement
    if condition is True:
        if firstRootElement is False:
            sys.stdout.write(',')
        firstRootElement = False
def conditionalPrint(condition, text):
    if condition is True:
        sys.stdout.write(text)


# Tools for handling region names
codeLinkedEventParser = re.compile(r'(/phylanx/[^\$]*)\$([^\$]*)\$([^\$]*)\$([^\$]*)\$([^\$]*)')

def addRegionChild(parent, child):
    if 'parent' in regions[child]:
        assert regions[child]['parent'] == parent
    else:
        regions[child]['parent'] = parent
    regions[parent]['children'].add(child)

def dumpRegion(regionName, parent=None):
    if regionName in regions:
        return
    regions[regionName] = { 'children': set() }
    if parent is not None:
        addRegionChild(parent, regionName)
        regions[regionName]['parent'] = parent
    codeLinkedEvent = codeLinkedEventParser.match(regionName)
    if codeLinkedEvent is not None:
        regions[regionName]['name'] = codeLinkedEvent[1]
        regions[regionName]['token1'] = codeLinkedEvent[2]
        regions[regionName]['token2'] = codeLinkedEvent[3]
        regions[regionName]['line'] = codeLinkedEvent[4]
        regions[regionName]['char'] = codeLinkedEvent[5]

# Tools for handling the tree
treeModeParser = re.compile(r'Tree information for function:')

regions = {}
def dumpTree(node, spaces, parent=None):
    # Create the hashed region object
    regionName = node.name.strip()
    assert regionName not in regions
    dumpRegion(regionName, parent)

    # Create the tree hierarchy
    conditionalPrint(args.tree, '\n' + spaces + '"name": "' + regionName + '",')
    if len(node.descendants) > 0:
        conditionalPrint(args.tree, '\n' + spaces + '"children": [')
        for index, child in enumerate(node.descendants):
            conditionalPrint(args.tree and index > 0, ',')
            conditionalPrint(args.tree, '\n' + spaces + '  {')
            dumpTree(child, spaces + '    ', parent=regionName)
            conditionalPrint(args.tree, '\n' + spaces + '  }')
        conditionalPrint(args.tree, '\n' + spaces + ']')
    else:
        conditionalPrint(args.tree, '\n' + spaces + '"children": []')

# Tools for handling the DOT graph
dotModeParser = re.compile(r'graph "[^"]*" {')
dotLineParser = re.compile(r'"([^"]*)" -- "([^"]*)";')

# Tools for handling the performance csv
perfModeParser = re.compile(r'primitive_instance,display_name,count,time,eval_direct')
perfLineParser = re.compile(r'"([^"]*)","([^"]*)",(\d+),(\d+),(-?1)')

# Parse stdout first (waits for it to finish before attempting to parse the OTF2 trace)
mode = None
for line in args.input:
    if mode is None:
        if treeModeParser.match(line):
            mode = 'tree'
            writeRootComma(args.tree)
            conditionalPrint(args.tree, '\n  "tree": {')
        elif dotModeParser.match(line):
            mode = 'dot'
        elif perfModeParser.match(line):
            mode = 'perf'
    elif mode == 'tree':
        dumpTree(newick.loads(line)[0], '    ')
        conditionalPrint(args.tree, '\n  }')
        mode = None
    elif mode == 'dot':
        dotLine = dotLineParser.match(line)
        if dotLine is not None:
            addRegionChild(dotLine[1], dotLine[2])
        else:
            mode = None
    elif mode == 'perf':
        perfLine = perfLineParser.match(line)
        if perfLine is not None:
            regionName = perfLine[1]
            dumpRegion(regionName, parent=None)     # TODO: assert regionName in regions fails sometimes... ask why
            regions[regionName]['display_name'] = perfLine[2]
            regions[regionName]['count'] = perfLine[3]
            regions[regionName]['time'] = perfLine[4]
            regions[regionName]['eval_direct'] = perfLine[5]
        else:
            mode = None
    else:
        # Should never reach this point
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

writeRootComma(args.events)
conditionalPrint(args.events, '\n  "events": {')
currentEvent = None
numEvents = 0
locations = {}

def dumpEvent(currentEvent, numEvents):
    if currentEvent is not None:
        numEvents += 1
        if numEvents % 10000 == 0:
            log('.', end=''),
        if numEvents % 100000 == 0:
            log('processed %i events' % numEvents)

        regionName = currentEvent['Region']
        dumpRegion(regionName, parent=None)     # TODO: figure out parent region based on GUIDs
        if 'eventCount' not in regions[regionName]:
            regions[regionName]['eventCount'] = 0
        regions[regionName]['eventCount'] += 1

        if args.ranges and (currentEvent['Event'] == 'ENTER' or currentEvent['Event'] == 'LEAVE'):
            # Add (and sort) the enter / leave event into its location list
            if not currentEvent['Location'] in locations:
                locations[currentEvent['Location']] = []
            heapq.heappush(locations[currentEvent['Location']], (currentEvent['Timestamp'], currentEvent))
        
        # Dump the event to the file
        conditionalPrint(args.events and numEvents > 1, ',')
        conditionalPrint(args.events, '\n    ' + json.dumps(currentEvent))
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
conditionalPrint(args.events, '\n  }')
log('finished processing %i events' % numEvents)

# Combine the sorted enter / leave events into ranges
if args.ranges is True:
    writeRootComma(args.ranges)
    sys.stdout.write('\n  "ranges":{')
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
                    log('.', end=''),
                if numRanges % 100000 == 0:
                    log('processed %i ranges' % numRanges)
                if numRanges > 0:
                    sys.stdout.write(',')
                sys.stdout.write('\n    ' + json.dumps(currentRange))
                lastEvent = None
        # Make sure there are no trailing ENTER events
        assert lastEvent is None
    # Finish the ranges dict
    sys.stdout.write('\n  }')
    log('finished processing %i ranges' % numRanges)

# Output the regions
writeRootComma(True)
sys.stdout.write('\n  "regions": {')
for rIndex, (regionName, region) in enumerate(regions.items()):
    region['children'] = list(region['children'])
    if rIndex > 0:
        sys.stdout.write(',')
    sys.stdout.write('\n    "' + regionName + '": {')
    for aIndex, (attrName, attr) in enumerate(region.items()):
        if aIndex > 0:
            sys.stdout.write(',')
        sys.stdout.write('\n      "' + attrName + '": ' + json.dumps(attr))
    sys.stdout.write('\n    }')
sys.stdout.write('\n  }')

# Output the graph links
if args.region_links:
    writeRootComma(True)
    sys.stdout.write('\n  "region links": {')
    for pIndex, (parent, value) in enumerate(regions.items()):
        for cIndex, child in enumerate(value['children']):
            if pIndex > 0 or cIndex > 0:
                sys.stdout.write(',')
            sys.stdout.write('\n    { "source": "' + parent + '", "target": "' + child + '" }')
    sys.stdout.write('\n  }')

# Finish the JSON output
sys.stdout.write('\n}')