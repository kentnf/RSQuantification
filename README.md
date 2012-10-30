1. Install

Perl modular: PerlIO::gzip, ....
Tophat
Bowtie
Samtools

2. Usage

build file index using bowtie
$ bowite2-build /path/plant_genome.fa /path/plant_genome

make a list for all cleaned fasta files 
$ ls *.fa > list

Then remove file suffix in list 

align the cleaned reads to genome
$ tophat_pipeline_v2.pl -i list -d /path/plant_genome -s SS -l firststrand -p 8 -m 1 -n 1 -x -a

get mapping rate for each samples (*make sure the cleaned reads locate at the same folder of bam files)
$ get_uniq_mapped_read.pl -i list -s SS

get the raw count of gene expression
$ get_exp_raw.pl -i list -s SS -a /path/plant_gene_position

normalization
$ get_exp_rpkm -e -x exp_sense_adjust /path/plant_gene_length
