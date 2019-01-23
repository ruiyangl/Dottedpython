# Dottedpython
dotplot maker in python for fasta sequence pairwise dot plot matrix analysis
Ruiyang Li, 2018
cridit also goes to Arvis Sulovari, Mark Chaisson, David Gordon and Peter Audano

#########################
#########################
### HAPPTYDOT PIPLINE ###
#########################
#########################


This pipeline produces haplotype resolved contigs and characterizes using multiple comparison dotplot matrix.
Required files:
	1..bam file of PacBio reads
	2.phased .vcf file
	3..bed file of the regions
There are 3 parts of this pipeline:
	1.generate haplotype resolved contigs from partitioned reads and haplotype specific SNVs
	2.map haplotype resolved contigs to human reference and extract sequence for the regions of interest
	3.create dotplots and sequence characterization 

@@@@ NOTE @@@@
all the required modules, if exist, should be loaded automatically, except before running the step GENERATE HAPLOTYPE RESOLVED CONTIGS, 
you need to manually load the required modules

#########################################
## GENERATE HAPLOTYPE RESOLVED CONTIGS ##
#########################################

++++++++++++++++++++++
++ REQUIRED MODULES ++
++++++++++++++++++++++

$module load miniconda/4.3.21

+++++++++++++++++++++++++++++
++ RUN 21K REGIONS ON CLINT++
+++++++++++++++++++++++++++++

in your working directory:

$~dgordon/phasedsv/scripts/copy_files.sh
	
make a .fofn file contains the path to the bam
	-----------------------------------------------------------------------------------------------------
	|"/net/eichler/vol2/home/dgordon/phasedsv/test_data/data/bams.fofn":								|
	|/gpfs/eichler27/projects/STR_VNTR/nobackups/phasing/phasedsv/clint/bams/output/kamilah_merged.bam  |
	-----------------------------------------------------------------------------------------------------

edit HG00733.conf
	1. to make DEST refer to your current directory (“run”)
	2. update BAMS to point to the fofn file
	3. update VCF to point to the phased .vcf file
	-------------------------------------------------------------------------------------------------------------
	|"HG00733.conf":																							|
	|REF=/net/eichler/vol2/eee_shared/assemblies/GRCh38/no_unloc/GRCh38.fasta									|
	|BAMS=/net/eichler/vol27/projects/ruiyang_projects/nobackups/vntr_project/phasedsv/runclint/clint.bams.fofn	|
	|VCF=/gpfs/eichler27/projects/STR_VNTR/nobackups/phasing/10X/analysis/clint/phased_variants.vcf.gz			|
	|SAMPLE=clint 																								|
	|DEST=/net/eichler/vol1/home/ruiyangl/nobackups/vntr_project/phasedsv/run21k								|
	-------------------------------------------------------------------------------------------------------------

edit asm.0.run.sh
	1. "region count" (1-21442 in the example)
	2. bed file points to the .bed file conteain the region coordinates
		This bed file needs to be in the "chr:start-end" formate
		-----------------------------------------------------------------------
		|"/net/eichler/vol2/home/dgordon/phasedsv/test_data/data/regions.txt:"|
		|chr1.23439975-23499975												  |
		|chr1.23459975-23519975												  |
		|chr1.23479975-23539975												  |
		|chr1.23499975-23559975												  |
		|chr1.23519975-23579975												  |
		|chr1.23539975-23599975												  |
		|chr1.23559975-23619975												  |
		|chr1.23579975-23639975												  |
		|chr1.23599975-23659975												  |
		|chr1.23619975-23679975												  |
		|chr1.23639975-23699975												  |
		|chr1.23659975-23719975												  |
		|chr1.23679975-23739975												  |
		-----------------------------------------------------------------------
	--------------------------------------------------------------------------------------------------------------------------------------------------------
	|"asm.0.run.sh":                         																											   |
	|#!/usr/bin/env bash																																   |
	|#$ -t 1-21442 -tc 50 -S /bin/bash -V  -e ./log -o ./log -l mfree=2G -l h_rt=01:00:00 -pe serial 4 -soft -l gpfsstate=0								   |
	|mkdir -p log																																		   |
	|source /net/eichler/vol2/home/dgordon/phasedsv/180702/phasedsv/config.sh 																			   |
	|/net/eichler/vol2/home/dgordon/phasedsv/180702/phasedsv/local_assembly/grid_scripts/../RunTiledAssembly.sh `awk 									   |
	|"NR == $SGE_TASK_ID" /net/eichler/vol27/projects/ruiyang_projects/nobackups/vntr_project/phasedsv/run21k/21444_regions_phaseSVformat.bed` HG00733.conf|
	--------------------------------------------------------------------------------------------------------------------------------------------------------
	
	@@@@@@ NOTE @@@@@@@
	There are only 21442 regions in the .bed file, however the name of the file inicates that there should be 214444 files


edit phasedsv.json
	1. "BAMs" to point to the same fofn file in HG00733.conf
	2. "sample" to point to the correct sample name

	--------------------------------------------------------------------------------------------------------------------------
	|"phasedsv.json":																										 |																				
	|{																														 |
	|                "ref": "/net/eichler/vol2/eee_shared/assemblies/GRCh38/no_unloc/GRCh38.fasta",							 |
	|          "sample": "Clint",																							 |
	|        "bams": "/net/eichler/vol27/projects/ruiyang_projects/nobackups/vntr_project/phasedsv/runclint/clint.bams.fofn",|
	|          "cov_cutoff": 3,																								 |
	|          "inversions" : "inversions/inversions.bed",																	 |
	|          "tr_cluster_size" : 6,																						 |
	|      "depth" : 40																										 |
	|}																														 |
	--------------------------------------------------------------------------------------------------------------------------

$qsub -cwd asm.0.run.sh

// To monitor the jobs
$watch -n 2 qstat -u <your user name>

after the jobs are done, the following files should be in <running directory>/asm/:
	1. h1_SV.fasta
	2. h2_SV.fasta
	3. assmebled_region.txt
	4. failed_reigon.txt

	@@@@ NOTE @@@@
	not all of the names in the failed_region.txt are valid regions, all the files that does not meet the requirments would appear to be filed region
	some intermediate files might also exist in this directory ignore them



####################################
## MAP HAPLOTIGS TO THE REFERACNE ##
####################################

++++++++++++++++++++++
++ REQUIRED MODULES ++
++++++++++++++++++++++

modules modules-init modules-gs/prod modules-eichler
mpc/0.8.2 mpfr/3.1.0 gmp/5.0.2 gcc/4.9.1 R/3.5.0
bedtools/2.25.0 bioawk/20110810 samtools/1.4
minimap2 miniconda/4.3.21


++++++++++++++++++++++
++ RUN MAP HAPLOTIG ++
++++++++++++++++++++++

In your working directory:

$/net/eichler/vol27/projects/ruiyang_projects/nobackups/vntr_project/phasedsv/map_haplotig/copy_map_haplotigs.sh

modify map_haplotig.json:
	1. change input 'target_region', this should be a .bed file in Tab Separated Values format
	2. running directory 'run_dir'
		the running directory should be the same running directory as the Phased-SV working directory

$bash run_map

when it is finished there should be 

###################################################
## CREATE DOTPLOTS AND SEQUENCE CHARACTERIZATION ##
###################################################

++++++++++++++++++++++
++ REQUIRED MODULES ++
++++++++++++++++++++++

$module load miniconda/4.3.21
$export DRMAA_LIBRARY_PATH=/opt/uge/lib/lx-amd64/libdrmaa.so.1.0
$module load trf/4.07b
$module load samtools/1.9

+++++++++++++++++++++++++++
++ RUN DOTPLOT GENERATOR ++
+++++++++++++++++++++++++++

in your working directory:

$~ruiyangl/nobackups/vntr_project/dotplots/copy_dotplot.sh

edit dotplots.json:
	1. "wind_size" to be the sliding window size for dotplots
	2. "output_dir"
	3. "tmp_dir" to point to where the intermediate files would be stored
	4. region file which should be the same one that is used in step one
	5. region count
	6. sample list, which should have the order that you want them to appear in the dotplot
		the first rule of the pipeline generates a fasta file of the sequence from the reference of the given regions in the reference (a lot of prepositions here sorry about that)
		it wouldn't be used in generating dotplots if you don't put it in the sample list here

	-----------------------------------------------------------------------------------------------------------------------------
	|"dotplots.json":																											|
	|{																															|
	|        "referance":"/net/eichler/vol2/eee_shared/assemblies/hg38/ucsc.hg38.fasta",										|
	|        "refgene":"/net/eichler/vol27/projects/ruiyang_projects/nobackups/vntr_project/dotplots/data/refgene/refgene.txt",	|
	|        "wind_size":10,																									|
	|        "output_dir":"results/644_human_specific",																			|
	|        "tmp_dir":"temp/644_human_specific",																				|
	|        "region_file":"data/644_human_specific/region.txt",																|
	|        "region_count":644,																								|
	|        "sample_list":[																									|
	|        "temp/644_human_specific/GRCh38.fasta",																			|
	|        "data/644_human_specific/HG00514.h0.fasta",																		|	
	|        "data/644_human_specific/HG00514.h1.fasta",																		|
	|        "data/644_human_specific/HG00733.h0.fasta",																		|
	|        "data/644_human_specific/HG00733.h1.fasta",																		|
	|        "data/644_human_specific/NA19240.h0.fasta",																		|
	|        "data/644_human_specific/NA19240.h1.fasta",																		|
	|        "data/644_human_specific/yoruban_hsa.fasta",																		|	
	|        "data/644_human_specific/chm13_hsa.fasta",																			|
	|        "data/644_human_specific/clint_ptr.fasta",																			|
	|        "data/644_human_specific/clint_ptr_h1.fasta",																		|	
	|        "data/644_human_specific/clint_ptr_h2.fasta",																		|	
	|        "data/644_human_specific/susie_ggo.fasta",																			|
	|        "data/644_human_specific/kamilah_ggo.h1.fasta",																	|			
	|        "data/644_human_specific/kamilah_ggo.h2.fasta",																	|			
	|        "data/644_human_specific/susie_pab.fasta",																			|	
	|        "data/644_human_specific/susie_pab_h1.fasta",																		|		
	|        "data/644_human_specific/susie_pab_h2.fasta"																		|		
	|        ]																													|
	|}																															|
	-----------------------------------------------------------------------------------------------------------------------------

Customize the colors and names for the dotplots of each species in dotplots/Input.py
by changing the list variable HIGHLIGHT

$bash run_dotplots_snake.sh

when it is finished the following files should be in the output directory:
	1. dotplot_<region>.pdf
	2. fasta_gene_list.csv
	3. motif_stats.tsv
