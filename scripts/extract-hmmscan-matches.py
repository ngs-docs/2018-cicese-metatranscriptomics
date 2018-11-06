#! /usr/bin/env python
import argparse, screed, sys

parser = argparse.ArgumentParser()
parser.add_argument('names')
parser.add_argument('seqfile')
args = parser.parse_args()

names = set([ x.strip() for x in open(args.names) ])
found = set()

for record in screed.open(args.seqfile):
    shortname = record.name.split()[0]
    if shortname in names:
        print('>{}\n{}'.format(record.name, record.sequence))
        found.add(shortname)


sys.stderr.write('found {} of {} ({} missing)\n'.format(len(found), len(names), len(names - found)))
