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
multiqc_container = 'docker://ewels/multiqc:v1.11'
fastqc_container = 'docker://biocontainers/fastqc:v0.11.9_cv8'

#########
# RULES #
#########

rule target:
    input:
     expand('output/mh_salmon/{sample}_quant/quant.sf', sample=all_samples),
     expand('output/mh_old_salmon/{sample}_quant/quant.sf', sample=all_samples),
     'output/multiqc/multiqc_report.html',
     'output/old_multiqc/multiqc_report.html',
     #expand('output/deseq2/tissue_LRT/{tissue}/{tissue}_GO_enrichment.pdf', tissue=["head", "abdo", "thorax", "venom", "ovary"]),
     #'output/blast/viral_genes/blastn.outfmt6',
     #'output/blast/crawford_venom/blastn.outfmt6'

#######################################
## blast seq.s against transcriptome ##
#######################################

##where are venom genes expressed?
rule blast_crawford_seq:
    input:
        crawford = 'data/crawford_seq/mh_venom_nt.fasta',
        db = 'output/blast/transcriptome_blastdb/mh_transcriptome.nhr'
    output:
        blast_res = 'output/blast/crawford_venom/blastn.outfmt6'
    params:
        db = 'output/blast/transcriptome_blastdb/mh_transcriptome'
    threads:
        20
    log:
        'output/logs/blast_crawford_seq.log'
    shell:
        'blastn '
        '-query {input.crawford} '
        '-db {params.db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blast_res} '
        '2>{log}' 

rule tissue_specific_GO_enrichment:
    input:
        tissue_degs_file = 'output/deseq2/tissue_LRT/{tissue}/{tissue}_sp_LRT_annots.csv'
    output:
        enrichment_table = 'output/deseq2/tissue_LRT/{tissue}/{tissue}_GO_enrichment.csv',
        GO_plot = 'output/deseq2/tissue_LRT/{tissue}/{tissue}_GO_enrichment.pdf'
    threads:
        10
    log:
        'output/logs/{tissue}_specific_GO_enrichment.log'
    script:
        'src/LRT/tissue_specific_GO_enrichment.R'

##############################
## map to mh transcriptome ##
##############################

rule old_multiqc:
    input:
        expand('output/mh_old_salmon/{sample}_quant/quant.sf', sample=all_samples)
    output:
        'output/old_multiqc/multiqc_report.html'
    params:
        outdir = 'output/old_multiqc',
        indirs = ['output/mh_old_salmon', 'output/fastqc']
    log:
        'output/logs/old_multiqc.log'
    container:
        multiqc_container
    shell:
        'multiqc -f ' ##force to write over old output
        '-o {params.outdir} '
        '{params.indirs} '
        '2> {log}'

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

#######################
## old transcriptome ##
#######################

rule mh_old_salmon_quant:
    input:
        index_output = 'output/mh_old_salmon/transcripts_index/refseq.bin',
        trimmed_r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
        trimmed_r2 = 'output/bbduk_trim/{sample}_r2.fq.gz'
    output:
        quant = 'output/mh_old_salmon/{sample}_quant/quant.sf'
    params:
        index_outdir = 'output/mh_old_salmon/transcripts_index',
        outdir = 'output/mh_old_salmon/{sample}_quant'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/mh_gg_salmon/salmon_logs/mh_old_salmon_quant_{sample}.log'
    shell:
        'salmon quant '
        '-i {params.index_outdir} '
        '-l ISR '
        '-1 {input.trimmed_r1} '
        '-2 {input.trimmed_r2} '
        '-o {params.outdir} '
        '-p {threads} '
        '&> {log}'

rule mh_old_salmon_index:
    input:
        transcriptome_length_filtered = 'data/diff_sample_mh_transcriptome/trinity_filtered_isoforms/isoforms_by_length.fasta'
    output:
        'output/mh_old_salmon/transcripts_index/refseq.bin'
    params:
        outdir = 'output/mh_old_salmon/transcripts_index'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/logs/mh_old_salmon_index.log'
    shell:
        'salmon index '
        '-t {input.transcriptome_length_filtered} '
        '-i {params.outdir} '
        '-p {threads} '
        '&> {log}'

###########################
## de novo transcriptome ##
###########################

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
