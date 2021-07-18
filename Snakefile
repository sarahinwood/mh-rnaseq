#!/usr/bin/env python3
import peppy

###########
# GLOBALS #
###########

##this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
##can now use this to generate list of all samples
all_samples = pep.sample_table['sample_name']

#containers
salmon_container = 'docker://combinelab/salmon:1.5.1'
kraken_container = 'docker://staphb/kraken2:2.1.2-no-db'
multiqc_container = 'docker://ewels/multiqc:v1.11'
fastqc_container = 'docker://biocontainers/fastqc:v0.11.9_cv8'

#########
# RULES #
#########

rule target:
    input:
     expand('output/mh_salmon/{sample}_quant/quant.sf', sample=all_samples),
     'output/multiqc/multiqc_report.html',
     #expand('output/kraken/kraken_{sample}_out.txt', sample=all_samples)

###############################################
## look at mapping to genome - STAR probably ##
###############################################

##############################
## map to mh transcriptome ##
##############################

rule multiqc:
    input:
        expand('output/mh_salmon/{sample}_quant/quant.sf', sample=all_samples)
    output:
        'output/multiqc/multiqc_report.html'
    params:
        outdir = 'output/multiqc',
        indirs = ['output/mh_salmon', 'output/fastqc']
    log:
        'output/logs/multiqc.log'
    container:
        multiqc_container
    shell:
        'multiqc -f ' ##force to write over old output
        '-o {params.outdir} '
        '{params.indirs} '
        '2> {log}'

rule mh_salmon_quant:
    input:
        index_output = 'output/mh_salmon/transcripts_index/refseq.bin',
        trimmed_r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
        trimmed_r2 = 'output/bbduk_trim/{sample}_r2.fq.gz'
    output:
        quant = 'output/mh_salmon/{sample}_quant/quant.sf'
    params:
        index_outdir = 'output/mh_salmon/transcripts_index',
        outdir = 'output/mh_salmon/{sample}_quant'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/mh_salmon/salmon_logs/mh_salmon_quant_{sample}.log'
    shell:
        'salmon quant '
        '-i {params.index_outdir} '
        '-l ISR '
        '-1 {input.trimmed_r1} '
        '-2 {input.trimmed_r2} '
        '-o {params.outdir} '
        '-p {threads} '
        '&> {log}'

rule mh_salmon_index:
    input:
        transcriptome_length_filtered = 'data/mh-transcriptome/output/trinity_filtered_isoforms/isoforms_by_length.fasta'
    output:
        'output/mh_salmon/transcripts_index/refseq.bin'
    params:
        outdir = 'output/mh_salmon/transcripts_index'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/logs/mh_salmon_index.log'
    shell:
        'salmon index '
        '-t {input.transcriptome_length_filtered} '
        '-i {params.outdir} '
        '-p {threads} '
        '&> {log}'

##kraken db from https://benlangmead.github.io/aws-indexes/k2 as build std was failing
rule kraken:
    input:
        r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
        r2 = 'output/bbduk_trim/{sample}_r2.fq.gz',
        db = 'bin/db/kraken_std'
    output:
        out = 'output/kraken/kraken_{sample}_out.txt',
        report = 'output/kraken/reports/kraken_{sample}_report.txt'
    log:
        'output/logs/kraken/kraken_{sample}.log'
    threads:
        20
    singularity:
        kraken_container
    shell:
        'kraken2 '
        '--threads {threads} '
        '--db {input.db} '
        '--paired '
        '--output {output.out} '
        '--report {output.report} '
        '--use-names '
        '{input.r1} {input.r2} '
        '&> {log}'