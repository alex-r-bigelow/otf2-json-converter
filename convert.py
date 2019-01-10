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

firstLineParser = re.compile('^(\S+)\s+(\d+)\s+(\d+)\s+Region: "([^"]*)" <(\d+)>$')
secondLineParser = re.compile('^\s+ADDITIONAL ATTRIBUTES: (.*)$')
attrSplitter = re.compile('\), \(')
attrParser = re.compile('\(?"([^"]*)" <\d+>; [^;]*; ([^\)]*)')

outfile.write('{')

currentEvent = None
numEvents = 0

def dumpEvent(currentEvent, numEvents):
    if currentEvent is not None:
        outfile.write('\n');
        outfile.write(json.dumps(currentEvent))
        numEvents += 1
        if numEvents % 100 == 0:
            print('.', end=''),
        if numEvents % 1000 == 0:
            print('processed %i events' % numEvents)
    return numEvents

for line in otfPrint.stdout:
    line = line.decode()
    firstLineMatch = firstLineParser.match(line)
    secondLineMatch = secondLineParser.match(line)
    if currentEvent is None and firstLineMatch is None:
        continue

    if firstLineMatch:
        numEvents = dumpEvent(currentEvent, numEvents)
        currentEvent = {}
        currentEvent['Event'] = firstLineMatch.group(1)
        currentEvent['Location'] = int(firstLineMatch.group(2))
        currentEvent['Timestamp'] = int(firstLineMatch.group(3))
        currentEvent['Region'] = firstLineMatch.group(4) + ' ' + firstLineMatch.group(5)
    else:
        assert currentEvent is not None and secondLineMatch is not None
        attrList = attrSplitter.split(secondLineMatch.group(1))
        for attrStr in attrList:
            attr = attrParser.match(attrStr)
            assert attr is not None
            currentEvent[attr.group(1)] = attr.group(2)
dumpEvent(currentEvent, numEvents)
outfile.write('\n}')

outfile.close()
