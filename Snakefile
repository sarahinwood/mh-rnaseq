#!/usr/bin/env python3
import peppy

###########
# GLOBALS #
###########

## sample list
# this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
# can now use this to generate list of all samples
all_samples = pep.sample_table['sample_name']

## tissue list
all_tissues = ["Head", "Thorax", "Abdomen", "Ovary", "Venom"]

#containers
salmon_container = 'docker://combinelab/salmon:1.5.1'
multiqc_container = 'docker://ewels/multiqc:v1.11'
fastqc_container = 'docker://biocontainers/fastqc:v0.11.9_cv8'
blast_container = 'docker://ncbi/blast:2.12.0'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
bioconductor_container = 'library://sinwood/bioconductor/bioconductor_3.12:0.0.1'


#########
# RULES # doesn't need to re-run meiosis sp
#########

rule target:
    input:
        ### QC ###
     #'output/02_multiqc/multiqc_report.html',
        ### DESeq2 analysis ###
     'output/03_deseq2/PCA/PCA_tissue.pdf',
     #expand('output/03_deseq2/tissue_itWT/{tissue}/{tissue}_sp_all_annots.csv', tissue=all_tissues),
     #'output/03_deseq2/tissue_LRT/sig_degs.csv',
     #expand('output/03_deseq2/tissue_itWT/{tissue}/{tissue}_enrich_plot.pdf', tissue=["Head", "Thorax", "Venom", "Ovary"]), ##abdo has no enriched PFAM or GO  terms so R script fails
        ### meiosis-sp ###
     #'output/03_deseq2/meiosis_sp/meiosis_heatmap.pdf',
        ### Blast searches ###
         # ID transcriptome hits for Crawford genes & blast those best hits
     #'output/04_blast/crawford_transcriptome/best_transcriptome_nr_blast.csv',
         # blast crawford genes as is
     #'output/04_blast/crawford_nr/transcripts_best_nrblast.csv',
        # blast venom DEGs with signalp
     #'output/04_blast/venom_signalp/venom_signalp_nr_blastx.outfmt6',
     'supp_tables/venom_blast/new_nr_blastx.outfmt6'

## blast venom genes for paper table



rule blastx_genes:
    input:
        unann_deg_transcripts = 'supp_tables/venom_blast/transcripts_new.fasta'
    output:
        blastx_res = 'supp_tables/venom_blast/new_nr_blastx.outfmt6'
    params:
        blast_db = 'bin/blast_db/nr/nr'
    threads:
        50
    singularity:
        blast_container
    log:
        'supp_tables/venom_blast/venom_blastx.log'
    shell:
        'blastx '
        '-query {input.unann_deg_transcripts} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blastx_res} '
        '2> {log}'

rule filter_genes:
    input:
        mh_transcriptome = 'data/mh-transcriptome/output/trinity_filtered_isoforms/isoforms_by_length.fasta',
        transcript_ids = 'supp_tables/venom_blast/new_to_blast.txt'
    output:
        transcripts = 'supp_tables/venom_blast/transcripts_new.fasta'
    singularity:
        bbduk_container
    log:
        'supp_tables/venom_blast/filter_genes.log'
    shell:
        'filterbyname.sh '
        'in={input.mh_transcriptome} '
        'include=t '
        'names={input.transcript_ids} '
        'substring=name '
        'out={output.transcripts} '
        '&> {log}'

###############################################
## blast venom degs with signalp for nr hits ##
###############################################

rule trinotate_signalp_venom_blastx:
    input:
        unann_deg_transcripts = 'output/04_blast/venom_trinotate_signalp/transcripts.fasta'
    output:
        blastx_res = 'output/04_blast/venom_signalp/venom_signalp_nr_blastx.outfmt6'
    params:
        blast_db = 'bin/blast_db/nr/nr'
    threads:
        50
    singularity:
        blast_container
    log:
        'output/logs/trinotate_signalp_venom_blastx.log'
    shell:
        'blastx '
        '-query {input.unann_deg_transcripts} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blastx_res} '
        '2> {log}'

rule filter_venom_signalp_fasta:
    input:
        mh_transcriptome = 'data/mh-transcriptome/output/trinity_filtered_isoforms/isoforms_by_length.fasta',
        transcript_ids = 'output/04_blast/venom_signalp/deg_ids_signalp.txt'
    output:
        transcripts = 'output/04_blast/venom_signalp/transcripts.fasta'
    singularity:
        bbduk_container
    log:
        'output/logs/filter_trinotate_signalp_venom.log'
    shell:
        'filterbyname.sh '
        'in={input.mh_transcriptome} '
        'include=t '
        'names={input.transcript_ids} '
        'substring=name '
        'out={output.transcripts} '
        '&> {log}'

# filter out any venom DEGs with signalP and blast all of them against nr
# 82 with signalp
# 66 of those no nr blast search yet
rule filter_venom_signalp:
    input:
        venom_file = 'output/03_deseq2/tissue_itWT/Venom/Venom_sp_all_annots.csv'
    output:
        deg_ids = 'output/04_blast/venom_signalp/deg_ids_signalp.txt'
    log:
        'output/logs/blast/filter_venom_signalp.log'
    singularity:
        bioconductor_container
    script:
        'src/blast_res/venom_sigP/filter_venom_signalp.R'

################################
## blast crawford transcripts ##
################################

#### blast best crawford transcripts against nr db #####

rule crawford_transcripts_blast_res:
    input:
        res_file = 'output/04_blast/crawford_transcriptome/transcript_nrblastn.outfmt6'
    output:
        all_res = 'output/04_blast/crawford_transcriptome/all_transcriptome_nr_blast.csv',
        best_hits = 'output/04_blast/crawford_transcriptome/best_transcriptome_nr_blast.csv'
    log:
        'output/logs/blast/crawford_transcripts_blast_res.log'
    singularity:
        bioconductor_container
    script:
        'src/blast_res/crawford/crawford_transcriptome_hits_blast_res.R'

rule nrblast_crawford_transcripts:
    input:
        transcripts = 'output/04_blast/crawford_transcriptome/transcripts.fasta',
        db = 'output/04_blast/transcriptome_blastdb/mh_transcriptome.nhr'
    output:
        blast_res = 'output/04_blast/crawford_transcriptome/transcript_nrblastn.outfmt6'
    params:
        blast_db = 'bin/blast_db/nr/nr'
    threads:
        20
    singularity:
        blast_container
    log:
        'output/logs/blast/nrblast_crawford_transcripts.log'
    shell:
        'blastx '
        '-query {input.transcripts} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blast_res} '
        '2>{log}'

rule filter_crawford_transcripts:
    input:
        mh_transcriptome = 'data/mh-transcriptome/output/trinity_filtered_isoforms/isoforms_by_length.fasta',
        transcript_ids = 'output/04_blast/crawford_transcriptome/crawford_transcript_ids.txt'
    output:
        transcripts = 'output/04_blast/crawford_transcriptome/transcripts.fasta'
    singularity:
        bbduk_container
    log:
        'output/logs/filter_crawford_transcripts.log'
    shell:
        'filterbyname.sh '
        'in={input.mh_transcriptome} '
        'include=t '
        'names={input.transcript_ids} '
        'substring=name '
        'out={output.transcripts} '
        '&> {log}'

###### blast crawford genes against transcriptome ######

rule crawford_transcriptome_res:
    input:
        res_file = 'output/04_blast/crawford_transcriptome/blastn.outfmt6',
        trinotate_file = 'data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv',
        venom_degs_file = 'output/03_deseq2/tissue_itWT/Venom/Venom_sp_LRT_annots.csv',
        transcript_lengths_file = 'data/mh-transcriptome/output/trinity_abundance/RSEM.isoforms.results',
        crawford_lengths_file = 'data/crawford_seq/crawford_gene_lengths.csv',
        salmon_tpm_file = 'output/03_deseq2/salmon_TPM.csv'
    output:
        all_res = 'output/04_blast/crawford_transcriptome/transcripts_all_nrblast.csv',
        best_res = 'output/04_blast/crawford_transcriptome/transcripts_best_blastn.csv',
        crawford_venom_degs = 'output/04_blast/crawford_transcriptome/crawford_venom_degs.csv',
        best_hit_transcript_ids = 'output/04_blast/crawford_transcriptome/crawford_transcript_ids.txt'
    log:
        'output/logs/blast/crawford_transcriptome_res.log'
    singularity:
        bioconductor_container
    script:
        'src/blast_res/crawford/crawford_transcriptome_res.R'

rule blast_crawford_seq_transcriptome:
    input:
        crawford = 'data/crawford_seq/mh_venom_nt.fasta',
        db = 'output/04_blast/transcriptome_blastdb/mh_transcriptome.nhr'
    output:
        blast_res = 'output/04_blast/crawford_transcriptome/blastn.outfmt6'
    params:
        blast_db = 'output/04_blast/transcriptome_blastdb/mh_transcriptome'
    threads:
        20
    singularity:
        blast_container
    log:
        'output/logs/blast/blast_crawford_seq_transcriptome.log'
    shell:
        'blastn '
        '-query {input.crawford} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blast_res} '
        '2>{log}'

rule mh_blast_db:
    input:
        transcriptome = 'data/mh-transcriptome/output/trinity_filtered_isoforms/isoforms_by_length.fasta'
    output:
        'output/04_blast/transcriptome_blastdb/mh_transcriptome.nhr'
    params:
        db = 'output/04_blast/transcriptome_blastdb/mh_transcriptome'
    threads:
        10
    log:
        'output/logs/blast/mh_blast_db.log'
    shell:
        'makeblastdb '
        '-in {input.transcriptome} '
        '-dbtype nucl '
        '-title mh_transcriptome '
        '-out {params.db} '
        '-parse_seqids '
        '2> {log}'

###### blast crawford genes as is ######

rule crawford_nrblast_res:
    input:
        res_file = 'output/04_blast/crawford_nr/blastx.outfmt6'
    output:
        all_res = 'output/04_blast/crawford_nr/transcripts_all_nrblast.csv',
        best_res = 'output/04_blast/crawford_nr/transcripts_best_nrblast.csv'
    log:
        'output/logs/blast/crawford_nrblast_res.log'
    singularity:
        bioconductor_container
    script:
        'src/blast_res/crawford/crawford_nrblast_res.R'

rule crawford_nrblast:
    input:
        'data/crawford_seq/mh_venom_nt.fasta'
    output:
        blastx_res = 'output/04_blast/crawford_nr/blastx.outfmt6'
    params:
        blast_db = 'bin/blast_db/nr/nr'
    threads:
        20
    singularity:
        blast_container
    log:
        'output/logs/blast/crawford_nrblast.log'
    shell:
        'blastx '
        '-query {input} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blastx_res} '
        '2> {log}'

######################
## meiosis sp genes ##
######################

rule meiosis_specific_gene_expression:
    input:
        ovary_LRT_file = 'output/03_deseq2/tissue_LRT/sig_degs.csv',
        ovary_LRT_dds_file = 'output/03_deseq2/tissue_LRT/all_mh_tissue_LRT.rds',
        meiosis_sp_genes_file = 'output/03_deseq2/meiosis_sp/meiosis_genes.csv',
        meiosis_sp_blast_res = 'output/03_deseq2/meiosis_sp/transcripts_best_nrblast.csv'
    output:
        heatmap = 'output/03_deseq2/meiosis_sp/meiosis_heatmap.pdf'
    singularity:
        bioconductor_container
    log:
        'output/logs/deseq2/meiosis_specific_genes.log'
    script:
        'src/meiosis_sp/meiosis_specific_gene_expression.R'

rule ID_meiosis_specific_genes:
    input:
        trinotate_file = 'data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv',
        ovary_itWT_file = 'output/03_deseq2/tissue_itWT/Ovary/Ovary_sp_annots.csv',
        dds_file = 'output/03_deseq2/tissue_LRT/all_mh_tissue_LRT.rds',
        blastx_res = 'output/03_deseq2/meiosis_sp/dmc1_rec8/dmc1_rec8_nr_blastx.outfmt6'
    output:
        meiosis_sp_genes = 'output/03_deseq2/meiosis_sp/meiosis_genes.csv',
        meiosis_ids = 'output/03_deseq2/meiosis_sp/meiosis_ids.txt'
    singularity:
        bioconductor_container
    log:
        'output/logs/deseq2/ID_meiosis_specific_genes.log'
    script:
        'src/meiosis_sp/meiosis_specific_gene_expression.R'

rule tissue_lrt:
    input:
        mh_dds_file = 'output/03_deseq2/mh_dds.rds'
    output:
        dds = 'output/03_deseq2/tissue_LRT/all_mh_tissue_LRT.rds',
        all_res = 'output/03_deseq2/tissue_LRT/all_res.csv',
        sig_res = 'output/03_deseq2/tissue_LRT/sig_degs.csv'
    singularity:
        bioconductor_container
    log:
        'output/logs/deseq2/tissue_lrt.log'
    script:
        'src/tissue_LRT.R'

rule meiosis_nrblast_res:
    input:
        res_file = 'output/03_deseq2/meiosis_sp/nr_blastx.outfmt6'
    output:
        all_res = 'output/03_deseq2/meiosis_sp/transcripts_all_nrblast.csv',
        best_res = 'output/03_deseq2/meiosis_sp/transcripts_best_nrblast.csv'
    log:
        'output/logs/blast/meiosis_nrblast_res.log'
    singularity:
        bioconductor_container
    script:
        'src/meiosis_sp/meiosis_nrblast_res.R'

rule blast_meiosis_genes:
    input:
        meiosis_transcripts = 'output/03_deseq2/meiosis_sp/meiosis_genes.fasta'
    output:
        blastx_res = 'output/03_deseq2/meiosis_sp/nr_blastx.outfmt6'
    params:
        blast_db = 'bin/blast_db/nr/nr'
    threads:
        50
    singularity:
        blast_container
    log:
        'output/logs/blast_meiosis_genes.log'
    shell:
        'blastx '
        '-query {input.meiosis_transcripts} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blastx_res} '
        '2> {log}'

rule filter_meiosis_genes:
    input:
        mh_transcriptome = 'data/mh-transcriptome/output/trinity_filtered_isoforms/isoforms_by_length.fasta',
        transcript_ids = 'output/03_deseq2/meiosis_sp/meiosis_ids.txt'
    output:
        meiosis_transcripts = 'output/03_deseq2/meiosis_sp/meiosis_genes.fasta'
    singularity:
        bbduk_container
    log:
        'output/logs/filter_meiosis_genes.log'
    shell:
        'filterbyname.sh '
        'in={input.mh_transcriptome} '
        'include=t '
        'names={input.transcript_ids} '
        'substring=name '
        'out={output.meiosis_transcripts} '
        '&> {log}'

# serch for DMC1 or REC8  - had no trinotate hits
rule nr_blast_dmc1_rec8_hits:
    input:
        meiosis_transcripts = 'output/03_deseq2/meiosis_sp/dmc1_rec8/dmc1_rec8_hits.fasta'
    output:
        blastx_res = 'output/03_deseq2/meiosis_sp/dmc1_rec8/dmc1_rec8_nr_blastx.outfmt6'
    params:
        blast_db = 'bin/blast_db/nr/nr'
    threads:
        50
    singularity:
        blast_container
    log:
        'output/logs/blast_meiosis_genes.log'
    shell:
        'blastx '
        '-query {input.meiosis_transcripts} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blastx_res} '
        '2> {log}'

rule filter_dmc1_rec8_hits:
    input:
        mh_transcriptome = 'data/mh-transcriptome/output/trinity_filtered_isoforms/isoforms_by_length.fasta',
        transcript_ids = 'output/03_deseq2/meiosis_sp/dmc1_rec8/dmc1_rec8_hit_ids.txt'
    output:
        meiosis_transcripts = 'output/03_deseq2/meiosis_sp/dmc1_rec8/dmc1_rec8_hits.fasta'
    singularity:
        bbduk_container
    log:
        'output/logs/filter_meiosis_genes.log'
    shell:
        'filterbyname.sh '
        'in={input.mh_transcriptome} '
        'include=t '
        'names={input.transcript_ids} '
        'substring=name '
        'out={output.meiosis_transcripts} '
        '&> {log}'

rule dmc1_rec8_hit_ids:
    input:
        res_file = 'output/03_deseq2/meiosis_sp/dmc1_rec8/dmc1_rec8_blastx.outfmt6'
    output:
        hit_ids = 'output/03_deseq2/meiosis_sp/dmc1_rec8/dmc1_rec8_hit_ids.txt'
    log:
        'output/logs/blast/dmc1_rec8_hit_ids.log'
    singularity:
        bioconductor_container
    script:
        'src/meiosis_sp/dmc1_rec8_hit_ids.R'

rule dmc1_rec8_transcriptome_blast:
    input:
        mh_transcriptome = 'data/mh-transcriptome/output/trinity_filtered_isoforms/isoforms_by_length.fasta',
        db = 'output/03_deseq2/meiosis_sp/dmc1_rec8/dmc1_rec8_db/dmc1_rec8.phr'
    output:
        blast_res = 'output/03_deseq2/meiosis_sp/dmc1_rec8/dmc1_rec8_blastx.outfmt6'
    params:
        db = 'output/03_deseq2/meiosis_sp/dmc1_rec8/dmc1_rec8_db/dmc1_rec8'
    threads:
        20
    singularity:
        blast_container
    log:
        'output/logs/blast/dmc1_rec8_transcriptome_blast.log'
    shell:
        'blastx '
        '-query {input.mh_transcriptome} '
        '-db {params.db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blast_res} '
        '2>{log}'

rule dmc1_rec8_blast_db:
    input:
        dmc1_rec8 = 'data/DMC1_REC8.fa'
    output:
        'output/03_deseq2/meiosis_sp/dmc1_rec8/dmc1_rec8_db/dmc1_rec8.phr'
    params:
        db = 'output/03_deseq2/meiosis_sp/dmc1_rec8/dmc1_rec8_db/dmc1_rec8'
    threads:
        10
    log:
        'output/logs/blast/mh_blast_db.log'
    shell:
        'makeblastdb '
        '-in {input.dmc1_rec8} '
        '-dbtype prot '
        '-title mh_transcriptome '
        '-out {params.db} '
        '-parse_seqids '
        '2> {log}'

#############################
## blast unann tissue DEGs ##
############################# 

rule merge_deg_annots_unann_blast:
    input:
        deg_file = 'output/03_deseq2/tissue_itWT/{tissue}/{tissue}_sp_annots.csv',
        blast_res_file = 'output/04_blast/unann_degs/blast_best_res.csv'
    output:
        degs_all_annots = 'output/03_deseq2/tissue_itWT/{tissue}/{tissue}_sp_all_annots.csv'
    singularity:
        bioconductor_container
    log:
        'output/logs/merge_deg_annots_unann_blast_{tissue}.log'
    script:
        'src/blast_res/unann/unann_merge_res_deg_lists.R'

rule unann_degs_blast_res:
    input:
        blast_res_file = 'output/04_blast/unann_degs/nr_blastx.outfmt6'
    output:
        blast_best_res = 'output/04_blast/unann_degs/blast_best_res.csv'
    singularity:
        bioconductor_container
    log:
        'output/logs/unann_des_blast_res.log'
    script:
        'src/blast_res/unann/unann_degs_blast_res.R'

rule unann_degs_blastx:
    input:
        unann_deg_transcripts = 'output/04_blast/unann_degs/unann_deg_transcripts.fasta'
    output:
        blastx_res = 'output/04_blast/unann_degs/nr_blastx.outfmt6'
    params:
        blast_db = 'bin/blast_db/nr/nr'
    threads:
        50
    singularity:
        blast_container
    log:
        'output/logs/unann_degs_blastx.log'
    shell:
        'blastx '
        '-query {input.unann_deg_transcripts} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blastx_res} '
        '2> {log}'

rule filter_unann_deg_transcripts:
    input:
        mh_transcriptome = 'data/mh-transcriptome/output/trinity_filtered_isoforms/isoforms_by_length.fasta',
        transcript_ids = 'output/04_blast/deg_ids_no_blastx.txt'
    output:
        unann_deg_transcripts = 'output/04_blast/unann_degs/unann_deg_transcripts.fasta'
    singularity:
        bbduk_container
    log:
        'output/logs/filter_unann_deg_transcripts.log'
    shell:
        'filterbyname.sh '
        'in={input.mh_transcriptome} '
        'include=t '
        'names={input.transcript_ids} '
        'substring=name '
        'out={output.unann_deg_transcripts} '
        '&> {log}'

rule id_degs_no_blastx:
    input:
        tissue_sp_annots = expand('output/03_deseq2/tissue_itWT/{tissue}/{tissue}_sp_annots.csv', tissue=all_tissues),
        sig_annots = 'output/03_deseq2/stage_WT/sig_annots.csv'
    output:
        deg_ids_no_blastx = 'output/04_blast/deg_ids_no_blastx.txt'
    singularity:
        bioconductor_container
    log:
        'output/logs/blast/id_degs_no_blastx.log'
    script:
        'src/blast_res/unann/id_degs_no_blastx.R'

##############
## stage DE ##
##############

rule stage_WT:
    input:
        mh_dds_file = 'output/03_deseq2/mh_dds.rds',
        trinotate_file = 'data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv',
        tissue_sp_annots = expand('output/03_deseq2/tissue_itWT/{tissue}/{tissue}_sp_annots.csv', tissue=all_tissues)
    output:
        stage_dds = 'output/03_deseq2/stage_WT/mh_stage_dds.rds',
        res_group = 'output/03_deseq2/stage_WT/res_group.csv',
        sig_annots = 'output/03_deseq2/stage_WT/sig_annots.csv',
        volcano = 'output/03_deseq2/stage_WT/volcano.pdf',
        venn_venom = 'output/03_deseq2/stage_WT/venn_venom.pdf',
        venn_ovary = 'output/03_deseq2/stage_WT/venn_ovary.pdf'
    singularity:
        bioconductor_container
    log:
        'output/logs/deseq2/stage_WT.log'
    script:
        'src/stage/stage_WT.R'

########################
## tissue specific DE ##
########################

rule plot_tissue_enrichment:
    input:
        pfam_file = 'output/03_deseq2/tissue_itWT/{tissue}/{tissue}_PFAM_enrichment.csv',
        go_file  = 'output/03_deseq2/tissue_itWT/{tissue}/{tissue}_GO_enrichment.csv'
    output:
        enrich_plot = 'output/03_deseq2/tissue_itWT/{tissue}/{tissue}_enrich_plot.pdf'
    singularity:
        bioconductor_container
    log:
        'output/logs/deseq2/plot_{tissue}_enrichment.log'
    script:
        'src/tissue_itWT/tissue_enrichment/plot_tissue_enrichment.R'

rule tissue_specific_term_enrichment:
    input:
        tissue_DEG_file = 'output/03_deseq2/tissue_itWT/{tissue}/{tissue}_sp_annots.csv',
        term_annot_table_file = 'output/03_deseq2/{term}_annots/{term}_annots.csv',
        term_to_gene_file = 'output/03_deseq2/{term}_annots/{term}_to_gene.csv',
        term_to_name_file = 'output/03_deseq2/{term}_annots/{term}_to_name.csv'
    output:
        enrichment_table = 'output/03_deseq2/tissue_itWT/{tissue}/{tissue}_{term}_enrichment.csv'
    singularity:
        bioconductor_container
    log:
        'output/logs/deseq2/{tissue}_specific_{term}_enrichment.log'
    script:
        'src/tissue_itWT/tissue_enrichment/tissue_specific_{wildcards.term}_enrichment.R'

rule term_to_gene:
    input:
        trinotate_file = 'data/mh-transcriptome/output/trinotate/trinotate/trinotate_annotation_report.txt'
    output:
        term_annot_table = 'output/03_deseq2/{term}_annots/{term}_annots.csv',
        term_to_gene = 'output/03_deseq2/{term}_annots/{term}_to_gene.csv',
        term_to_name = 'output/03_deseq2/{term}_annots/{term}_to_name.csv'
    singularity:
        bioconductor_container
    log:
        'output/logs/deseq2/{term}_term_to_gene.log'
    script:
        'src/tissue_itWT/tissue_enrichment/{wildcards.term}_to_geneID.R'

rule tissue_itWT_tests:
    input:
        mh_itWT_dds_file = 'output/03_deseq2/tissue_itWT/mh_itWT.rds',
        trinotate_file = 'data/mh-transcriptome/output/trinotate/sorted/longest_isoform_annots.csv',
        salmon_tpm_file = 'output/03_deseq2/salmon_TPM.csv'
    output:
        tissue_annots = 'output/03_deseq2/tissue_itWT/{tissue}/{tissue}_sp_annots.csv',
        itWT_venn = 'output/03_deseq2/tissue_itWT/{tissue}/{tissue}_itWT_venn.pdf',
        heatmap = 'output/03_deseq2/tissue_itWT/{tissue}/{tissue}_heatmap.pdf',
        clustered_heatmap = 'output/03_deseq2/tissue_itWT/{tissue}/{tissue}_clustered_heatmap.pdf'
    singularity:
        bioconductor_container
    log:
        'output/logs/deseq2/{tissue}_itWT_dds.log'
    script:
        'src/tissue_itWT/{wildcards.tissue}_itWT.R'

rule tissue_itWT_dds:
    input:
        mh_dds_file = 'output/03_deseq2/mh_dds.rds'
    output:
        mh_itWT_dds = 'output/03_deseq2/tissue_itWT/mh_itWT.rds'
    singularity:
        bioconductor_container
    log:
        'output/logs/deseq2/tissue_itWT_dds.log'
    script:
        'src/tissue_itWT/tissue_itWT_dds.R'

########################
## initial DGE set up ##
########################

rule PCA_analysis:
    input:
        mh_dds_file = 'output/03_deseq2/mh_dds.rds',
        sample_data_file = 'data/sample_table.csv'
    output:
        PCA_tissue = 'output/03_deseq2/PCA/PCA_tissue.pdf',
        PCA_batch = 'output/03_deseq2/PCA/PCA_replicate.pdf',
        PCA_seqrun = 'output/03_deseq2/PCA/PCA_seqrun.pdf'
    singularity:
        bioconductor_container
    log:
        'output/logs/deseq2/PCA_analysis.log'
    script:
        'src/QC/PCA_analysis.R'

rule make_mh_dds:
    input:
        gene2tx_file = 'data/mh-transcriptome/output/trinity/Trinity.fasta.gene_trans_map',
        salmon_output = expand('output/01_mh_salmon/{sample}_quant/quant.sf', sample=all_samples),
        sample_data_file = 'data/sample_table.csv'
    output:
        salmon_tpm = 'output/03_deseq2/salmon_TPM.csv',
        mh_dds = 'output/03_deseq2/mh_dds.rds',
    log:
        'output/logs/deseq2/make_dds.log'
    singularity:
        bioconductor_container
    script:
        'src/make_mh_dds.R'

#############################
## map to mh transcriptome ##
#############################

rule multiqc:
    input:
        expand('output/01_mh_salmon/{sample}_quant/quant.sf', sample=all_samples)
    output:
        'output/02_multiqc/multiqc_report.html'
    params:
        outdir = 'output/02_multiqc',
        indirs = ['output/01_mh_salmon', 'output/fastqc']
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
        index_output = 'output/01_mh_salmon/transcripts_index/refseq.bin',
        trimmed_r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
        trimmed_r2 = 'output/bbduk_trim/{sample}_r2.fq.gz'
    output:
        quant = 'output/01_mh_salmon/{sample}_quant/quant.sf'
    params:
        index_outdir = 'output/01_mh_salmon/transcripts_index',
        outdir = 'output/01_mh_salmon/{sample}_quant'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/01_mh_salmon/salmon_logs/mh_salmon_quant_{sample}.log'
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
        'output/01_mh_salmon/transcripts_index/refseq.bin'
    params:
        outdir = 'output/01_mh_salmon/transcripts_index'
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