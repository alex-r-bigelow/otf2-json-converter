#!/usr/bin/env python3
import sys
import argparse
import subprocess
import re
import json
import heapq

parser = argparse.ArgumentParser(description='Convert an OTF2 file to JSON')
parser.add_argument('-i', '--input', dest='input', default='test_data/OTF2_archive/APEX.otf2',
                    help='The input file (default=test_data/OTF2_archive/APEX.otf2)')
parser.add_argument('-e', '--events', dest='events',
                    help='A JSON file containing discrete events (e.g. ENTER and LEAVE are separate)')
parser.add_argument('-r', '--ranges', dest='ranges',
                    help='A JSON file containing combined ENTER and LEAVE events as single entries')

args = parser.parse_args()
if args.events is None and args.ranges is None:
    raise Exception('At least --events or --ranges is required')

otfPrint = subprocess.Popen(['otf2-print', args.input], stdout=subprocess.PIPE)
if args.events is not None:
    eventsFile = open(args.events, 'w')
    eventsFile.write('{')

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

currentEvent = None
numEvents = 0

locations = {}

def dumpEvent(currentEvent, numEvents):
    if currentEvent is not None:
        numEvents += 1
        if numEvents % 100 == 0:
            print('.', end=''),
        if numEvents % 1000 == 0:
            print('processed %i events' % numEvents)
        if args.events is not None:
            if numEvents > 1:
                eventsFile.write(',')
            eventsFile.write('\n  ' + json.dumps(currentEvent))
        if args.ranges is not None and (currentEvent['Event'] == 'ENTER' or currentEvent['Event'] == 'LEAVE'):
            if not currentEvent['Location'] in locations:
                locations[currentEvent['Location']] = []
            heapq.heappush(locations[currentEvent['Location']], (currentEvent['Timestamp'], currentEvent))
    return numEvents

for line in otfPrint.stdout:
    line = line.decode()
    eventLineMatch = eventLineParser.match(line)
    addAttrLineMatch = addAttrLineParser.match(line)
    if currentEvent is None and eventLineMatch is None:
        continue

    if eventLineMatch:
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
        assert currentEvent is not None and addAttrLineMatch is not None
        attrList = addAttrSplitter.split(addAttrLineMatch.group(1))
        for attrStr in attrList:
            attr = addAttrParser.match(attrStr)
            assert attr is not None
            currentEvent[attr.group(1)] = attr.group(2)
dumpEvent(currentEvent, numEvents)
if args.events is not None:
    eventsFile.write('\n}')
    eventsFile.close()

if args.ranges is not None:
    rangesFile = open(args.ranges, 'w')
    rangesFile.write('{')

    numRanges = 0
    for eventList in locations.values():
        lastEvent = None
        for timestamp, event in eventList:
            if event['Event'] == 'ENTER':
                assert lastEvent is None
                lastEvent = event
            elif event['Event'] == 'LEAVE':
                assert lastEvent is not None
                currentRange = {}
                for attr, value in event.items():
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

                numRanges += 1
                if numRanges % 100 == 0:
                    print('.', end=''),
                if numRanges % 1000 == 0:
                    print('processed %i ranges' % numRanges)
                if numRanges > 0:
                    rangesFile.write(',')
                rangesFile.write('\n  ' + json.dumps(currentRange))
                lastEvent = None
        assert lastEvent is None
    
    rangesFile.write('\n}')
    rangesFile.close()
