#!/usr/bin/env python3
import sys
import argparse
import subprocess
import re
import json

parser = argparse.ArgumentParser(description='Convert an OTF2 file to JSON')
parser.add_argument('-i', '--input', dest='input', default='test/als_csv_instrumented_MovieLens/APEX.otf2',
                    help='The input file (default=test/als_csv_instrumented_MovieLens/APEX.otf2)')
parser.add_argument('-o', '--output', dest='output', default='output.json',
                    help='The output file')

args = parser.parse_args()

otfPrint = subprocess.Popen(['otf2-print', args.input], stdout=subprocess.PIPE)
outfile = open(args.output, 'w')

eventLineParser = re.compile('^(\S+)\s+(\d+)\s+(\d+)\s+(.*)$')
attrParsers = {
  'ENTER': '(Region): "([^"]*)"',
  'LEAVE': '(Region): "([^"]*)"',
  'METRIC': 'Value: \("([^"]*)" <\d+>; [^;]*; ([^\)]*)',
  'MPI_SEND': '([^:]*): ([^,]*)',
  'MPI_RECV': '([^:]*): ([^,]*)'
}
addAttrLineParser = re.compile('^\s+ADDITIONAL ATTRIBUTES: (.*)$')
addAttrSplitter = re.compile('\), \(')
addAttrParser = re.compile('\(?"([^"]*)" <\d+>; [^;]*; ([^\)]*)')

outfile.write('{')

currentEvent = None
numEvents = 0

def dumpEvent(currentEvent, numEvents):
    if currentEvent is not None:
        outfile.write('\n  ' + json.dumps(currentEvent));
        numEvents += 1
        if numEvents % 100 == 0:
            print('.', end=''),
        if numEvents % 1000 == 0:
            print('processed %i events' % numEvents)
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
outfile.write('\n}')

outfile.close()
