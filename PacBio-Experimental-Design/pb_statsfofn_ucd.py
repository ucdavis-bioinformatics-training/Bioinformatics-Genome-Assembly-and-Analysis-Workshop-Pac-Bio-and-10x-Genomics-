#!/usr/bin/env python

'''
pb_statsfofn_ucd.py

generate read stats from pacbio fofn file

fofn format input, one line per fasta file
/share/dnat/rs2/161028_448/D02_1/Analysis_Results/m161031_221844_42145_c101058752550000001823247401061784_s1_p0.2.subreads.fasta
/share/dnat/rs2/161028_448/D02_1/Analysis_Results/m161031_221844_42145_c101058752550000001823247401061784_s1_p0.3.subreads.fasta

stats output file is in json format

'''
import sys
import os
import glob
import time
import re
import json
from collections import OrderedDict
from optparse import OptionParser  # http://docs.python.org/library/optparse.html

version = "1.0"
usage = "usage: %prog [options] -o output_filename fofn_file"
epilog = "The fofn_file can be provided as stdin. \
Specifying -o stdout can be used to put the ou tput to stdout."

parser = OptionParser(usage=usage, version="%prog " + str(version), epilog=epilog)
parser.add_option('-o', '--output', help="output filename, stdout is acceptable [default: %default]",
                  action="store", type="str", dest="output", default="input.fofn.stats")
parser.add_option('-d', '--id', help="ID for the stats output [default: %default]",
                  action="store", type="str", dest="id", default=None)

(options, args) = parser.parse_args()


def remove_comment(str1):
    loc = str1.find("#")  # identify if and where a # occurs
    if loc == -1:
        return str1  # if # not found return str1 as is
    str1 = str1[0:loc]  # trim of comment
    return str1.rstrip()  # remove any trailing whitespace and return result


def parse_pacbio_readid_rs2(str1):
    # >m160924_114322_42145_c101087062550000001823245703091787_s1_p0/8/0_35815 RQ=0.849
    rid, rq = str1.strip().split(' RQ=')  # strip out read quality
    rid, zmw, rlen = rid.split("/")
    rlen = rlen.split("_")
    rstart = int(rlen[0])
    rend = int(rlen[1])
    rlen = rend-rstart
    return (rid[1:], int(zmw), rstart, rend, rlen, float(rq))


if len(args) > 1:
    sys.stderr.write("program requires at most 1 argument\n")
    sys.stderr.write(usage + '\n')
    sys.exit()
elif len(args) == 1:
    infile = args[0]
    # Start opening input/output files:
    if not os.path.exists(infile):
        sys.stderr.write("Error, can't find input file %s\n" % infile)
        sys.exit()
    infofn = open(infile, 'r')
else:
    # reading from stdin
    infofn = sys.stdin

output = options.output

if output == "stdout":
    out = sys.stdout
else:
    out = open(output, 'w')

sid = options.id

fasta_files = 0

all_data = []

for line in infofn:
    # each line in the fofn file is one fasta file
    line = remove_comment(line)
    line = line.strip()
    if len(line) == 0:
        continue
    if not os.path.exists(line):
        sys.stderr.write("Error, can't find fasta file %s\n" % line)
        sys.exit()

    # process a fasta file
    fasta = open(line, 'r')
    sys.stderr.write("Processing file: %s\n" % line)
    # parse out run and cell position and cell barcode from file name
    pf = line.split("/Analysis_Results/")
    if len(pf) == 2:
        # agrees with UCD DNA Tech Core file path expectations
        pr = pf[0].split('/')
        run = pr[-2]
        m = re.match("m\d+_\d+_\d+_(c\d+)_s1_p0.\d.subreads.fasta",pf[1])
        if m:
            cell_position = pr[-1]
            cell_barcode = m.group(1)
        else:
            sys.stderr.write("Error, can't identify cell barcode in either filename, or first read\n")
            sys.exit()
    elif len(pf) == 1:
        run = None
        # try and extract cell_id from filename
        m = re.match(".+/m\d+_\d+_\d+_(c\d+)_s1_p0.\d.subreads.fasta",pf[0])
        if m:
            cell_position = None
            cell_barcode = m.group(1)
        else:
            # get cell id from first read id
            l1 = fasta.readline()
            m = re.match(">m\d+_\d+_\d+_(c\d+)_s1_p0/\d+/\d+_\d+ RQ=\d.\d+",l1)
            if m:
                cell_position = None
                cell_barcode = m.group(1)
            else:
                sys.stderr.write("Error, can't identify cell barcode in either filename, or first read\n")
                sys.exit()
            fasta.seek(0) # seek back to beginning of file

    else:
        sys.stderr.write("Error, can't identify cell barcode in either filename, or first read\n")
        sys.exit()            

    rcount = 0 # read count
    zmw = []
    rstart = []
    rend = []
    rlen = []
    rq =[]
    try:
        for read in fasta:
            if read[0] == '>':
                rstats = parse_pacbio_readid_rs2(read)
                zmw.append(rstats[1])
                rstart.append(rstats[2])
                rend.append(rstats[3])
                rlen.append(rstats[4])
                rq.append(rstats[5])
                rcount += 1
            else:
                continue
    finally:
        fasta.close()


    file_data = OrderedDict([("run_id", run), 
                             ("cell_position", cell_position),
                             ("cell_barcode", cell_barcode),
                             ("filename", line),
                             ("read_count", rcount),
                             ("zmw", zmw),
                             ("read_starts", rstart),
                             ("read_ends", rend),
                             ("read_lengths", rlen),
                             ("read_qualiies", rq)])

    all_data.append(file_data)
    file_data = None
    fasta_files += 1

file_keys = ["run_id","cell_position","cell_barcode","filename","read_count","zmw","read_starts","read_ends","read_lengths","read_qualities"]

stats_json = OrderedDict([
    ("id", sid),
    ("format", "UCD Pac Bio Fasta stats %s" % version),
    ("format_url", "https://github.com/ucdavis-bioinformatics/PacBio_GenomeStats"),
    ("generated_by", "pb_statsfofn_ucd.py"),
    ("date", time.strftime("%Y-%m-%dT%H:%M:%S")),
    ("type", "ucd pb stats"),
    ("source", "fasta"),
    ("number_of_files", fasta_files),
    ("file_keys", file_keys),
    ("file_data", all_data)])

sys.stderr.write("Writing JSON output to: %s\n" % output )
json.dump(stats_json, out)

infofn.close()
out.close()
