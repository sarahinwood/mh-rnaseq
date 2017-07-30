#!/usr/bin/env bash
set -eu
export LD_LIBRARY_PATH="bin/salmon/lib"
salmon="bin/salmon/bin/salmon"

"${salmon}" quant -i output/salmon/transcripts_index -l ISR -1 data/abdo_r1.fq.gz -2 data/abdo_r2.fq.gz -o output/salmon/abdo_quant -p 50

"${salmon}" quant -i output/salmon/transcripts_index -l ISR -1 data/head_r1.fq.gz -2 data/head_r2.fq.gz -o output/salmon/thorax_quant -p 50

"${salmon}" quant -i output/salmon/transcripts_index -l ISR -1 data/sting_r1.fq.gz -2 data/sting_r2.fq.gz -o output/salmon/thorax_quant -p 50