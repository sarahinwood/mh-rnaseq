#!/usr/bin/env python3
import pathlib2
import os
import pandas

#############
# FUNCTIONS #
#############

def resolve_path(x):
	return(str(pathlib2.Path(x).resolve(strict=False)))

def find_read_files(read_dir):
#Make list of files
	path_generator = os.walk(read_dir, followlinks = True)
	my_files = list((dirpath, filenames)
		for (dirpath, dirname, filenames)
		in path_generator)
#Make new dictionary & populate with files (flowcell = key)
	my_fastq_files = {}
	for dirpath, filenames in my_files:
		for filename in filenames:
			if filename.endswith('.fastq.gz'):
				my_flowcell = pathlib2.Path(dirpath).name
				my_fastq = str(pathlib2.Path(dirpath,filename))
				if my_flowcell in my_fastq_files:
					my_fastq_files[my_flowcell].append(my_fastq)
				else:
					my_fastq_files[my_flowcell]= []
					my_fastq_files[my_flowcell].append(my_fastq)
	return(my_fastq_files)

def sample_name_to_fastq(wildcards):
	sample_row = sample_key[sample_key['Sample_name'] == wildcards.sample]
	sample_id = sample_row.iloc[-1]['OGF_sample_ID']
	sample_flowcell = sample_row.iloc[-1]['Flow_cell']
	sample_all_fastq = [x for x in all_fastq[sample_flowcell]
						if '-{}-'.format(sample_id) in x]
	sample_r1 = sorted(list(x for x in sample_all_fastq
							if '_R1_' in os.path.basename(x)))
	sample_r2 = sorted(list(x for x in sample_all_fastq
							if '_R2_' in os.path.basename(x)))
	return({'r1': sample_r1, 'r2': sample_r2})

###########
# GLOBALS #
###########

read_dir = 'data/bbduk_trim'

sample_key_file = 'data/sample_key.csv'

#containers
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
salmon_container = 'shub://TomHarrop/singularity-containers:salmon_0.11.1'
star_container = 'shub://TomHarrop/singularity-containers:star_2.7.0c'

#########
# SETUP #
#########
# generate name to filename dictionary
all_fastq = find_read_files(read_dir)

sample_key = pandas.read_csv(sample_key_file)

all_samples = sorted(set(sample_key['Sample_name']))

#########
# RULES #
#########

rule target:
	input:
	 expand('output/2019_assembly_output/mh_salmon/{sample}_quant/quant.sf', sample = all_samples)

rule asw_salmon_quant:
	input:
		index_output = 'output/2019_assembly_output/mh_salmon/transcripts_index/hash.bin',
		left = 'data/bbduk_trim/{sample}_r1.fq.gz',
		right = 'data/bbduk_trim/{sample}_r2.fq.gz'
	output:
		quant = 'output/2019_assembly_output/mh_salmon/{sample}_quant/quant.sf'
	params:
		index_outdir = 'output/2019_assembly_output/mh_salmon/transcripts_index/',
		outdir = 'output/2019_assembly_output/mh_salmon/{sample}_quant'
	threads:
		20
	singularity:
		salmon_container
	log:
		'output/logs/salmon/mh_salmon_quant_2019_{sample}.log'
	shell:
		'salmon quant '
		'-i {params.index_outdir} '
		'-l ISR '
		'-1 {input.left} '
		'-2 {input.right} '
		'-o {params.outdir} '
		'-p {threads} '
		'&> {log}'

rule asw_salmon_index:
	input:
		transcriptome_length_filtered = 'data/isoforms_by_length.fasta'
	output:
		'output/2019_assembly_output/mh_salmon/transcripts_index/hash.bin'
	params:
		outdir = 'output/2019_assembly_output/mh_salmon/transcripts_index'
	threads:
		20
	singularity:
		salmon_container
	log:
		'output/logs/mh_salmon_index_2019.log'
	shell:
		'salmon index '
		'-t {input.transcriptome_length_filtered} '
		'-i {params.outdir} '
		'-p {threads} '
		'&> {log}'
