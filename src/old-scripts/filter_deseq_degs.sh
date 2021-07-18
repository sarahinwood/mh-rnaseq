#!/usr/bin/env bash

set -eu

abdohead_names="output/deseq2/abdovhead_sigresults.txt"
stinghead_names="output/deseq2/stingvhead_sigresults.txt"

bin/bbmap/filterbyname.sh \
in=data/isoforms_by_length.fasta \
include=t \
names="${abdohead_names}" \
out=output/deseq_filtered_degs/abdovhead_degs.fasta

bin/bbmap/filterbyname.sh \
in=data/isoforms_by_length.fasta \
include=t \
names="${stinghead_names}" \
out=output/deseq_filtered_degs/stingvhead_degs.fasta