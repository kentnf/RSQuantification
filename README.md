=================
1. Install
=================

Perl modular: PerlIO::gzip, ....\n
Tophat2
Bowtie2
samtools

=================
2. Usage
=================

2.1 aligning reads to reference 

build file index using bowtie
    $bowite2-build /path/plant_genome.fa /path/plant_genome

make a list for all cleaned fasta files 
    $ls *.fa > list

Then remove file suffix in list using VIM

Perform the fellowing analysis:
1) align the cleaned reads to genome
2) get No. of uniq mapped reads (make sure the cleaned reads locate at the same folder of bam files)
3) count the raw number for each gene
4) normalization using RPKM
    
    $tophat_pipeline_v2.pl -i list -d /path/plant_genome -s SS -l fr-firststrand -p 8 -m 1 -n 1 / 
                           -x gene_position_file -a gene_length

Perform the below analysis step by step:
1) align the cleaned reads to genome

    $tophat_pipeline_v2.pl -i list -d /path/plant_genome -s SS -l fr-firststrand -p 8 -m 1 -n 1

    Here is the usage of tophat_pipeline.pl
     -i  [String]   Name of the input list of read files (required)
     -s  [String]   Sequencing method for read files: (required)
            PE (paired-end);
            SE (single-end);
            PS (paired strand-specific);
            SS (single strand-specific);
     -d  [String]   Indexed database of genome sequences (required)
     -l  [String]   Library type, fr-unstrand, fr-firststrand, fr-secondstrand (required)
     -p  [Integer]  number of CPUs used for megablast clustering (default = 1)
     -r  [Integer]  mate-inner-dist (default = 100)
     -v  [Integer]  mate-std-dev (default = 20)
     -m  [Integer]  segment mismatches (default = 0)
     -n  [Integer]  initial read mismatches (default = 0)
     -e  [Integer]  setment length (default : half of read length; 17~25)
     -c  [Integer]  cycle (default = 0, range: 0,1,2)
     -a     (str)  gene and position file (optional)
     -x     (str)  gene and length file (optional)

    The required input includes list of read files, sequencing method, indexed genome
    sequences, and library type. The list is the prefix of read file. For example, if
    the read file name is WT_17D_rep2.fasta.gz, the list should be WT_17D_rep2. For
    paired-end read, the file name should be end with _1 or _2 before format info. For
    example, if the read file is test1_1.fq.gz and test1_2.fq.gz, the list is test1.
    Now, it supports gzip file with fasta or fastq format:
    ".gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz", ".fq", ".fq.gz", ".fastq", ".fastq.gz"

    The file name dose not affect the real format of read file, it a fastq file named
    as *.fasta, tophat_pipline.pl will identify the format automatically.

    The segment mismatches should be not shorter than initial read mismatches. If not,
    the tophat_pipline.pl will fix it.

    Please do not provide setment length, tophat_pipline.pl will set it automatically.

    If gene position and gene length files are provided, tophat_pipline.pl will generate
    raw count and RPKM info. Otherwise, we can generate these info follow the steps
    from 2 to 4.

2) get No. of uniq mapped reads (make sure the cleaned reads locate at the same folder of bam files)

    $get_uniq_mapped_read.pl -i list -s SS

    Here is the usage of get_uniq_mapped_num.pl.
	-i|list                 (str)  input read list file (required)
	-o|output               (str)  out file name   (default = uniq_mapped_num)
  	-s|sequencing-method    (str)  (required) (default = SS)
        	PE (paired-end);
        	SE (single-end);
        	PS (paired strand-specific);
        	SS (single strand-specific);

    The input list file is the same file of step 1.
    The sequencing-method must be same as step 1.

3) count the raw number for each gene

    $get_exp_raw.pl -i list -s SS -a gene_position_file

    Here is the usage of get_exp_raw.pl.
  	-i|list                (str)   input read list file (required)
  	-a|gene-position       (str)   gene position file (required)
  	-o|output              (str)   prefix of out files (default = exp; exp_sense_raw, exp_antisense_raw)
  	-s|sequencing-method   (str)   (required) (default = SS)
        	PE (paired-end);
        	SE (single-end);
        	PS (paired strand-specific);
        	SS (single strand-specific);
  	-h|?|help                help info

    The input list file is the same file of step 1.
    The sequencing-method must be same as step 1.

    The gene position file is a tab-delimited file, format is like below:
    chr  start  end  gene  strand  (the titile is not included)
    chr1 100    1000 Gene1 +

4) normalization using RPKM

    $get_exp_rpkm.pl -e exp_sense_raw -u uniq_mapped_read_num -x tomato_gene_length

    Here is the usage of get_exp_rpkm.pl.
  	-x|gene-length (str)   gene length
  	-e|exp-adjust  (str)   adjusted expression
  	-o|output      (str)   output (default = exp_rpkm)

    File uniq_mapped_num is the output of get_uniq_mapped_read.

    Gene length file is a tab-delimited file, format is like below:
    GeneID  Length (the titile is not included)
    Gene001 1000

5) Get correlation for each samples

    $corre.pl exp_rpkm > correlation.txt

2.2 Perform statistics analysis for differentially expressed genes

Prepare comparison file for pairwise comparison
    
    Samples: C_rep1, C_rep2, C_rep3, T_rep1, T_rep2, T_rep3
    Comparison_list File: C\tT\n
    
Pairwise statistics analysis using DESeq

    $DESeq_pipeline.pl  raw_count  rpkm_file  comparison_list  output

Pairwise statistics analysis using edgeR

    $edgeR_pipeline.pl  raw_count  rpkm_file  comparison_list  output

Prepare comparison file for time-series comparison

    Samples: T1h_rep1, T1h_rep2, T6h_rep1, T6h_rep2, T12h_rep1, T12h_rep2
    Comparison_list File: T1h\tT6h\tT12h\n

Time-series statistics analysis using DESeq and LIMMA

    $TSlimma_pipeline.pl  raw_count  rpkm_file  comparison_list  output

=================
3. Other Tools
=================

3.1 Prepare gene position file using GFF or GTF file
    In fact, using the representative transcriptome is better, 

3.2 Get unmapped reads

    $get_unmapped_read.pl -i list -s SS

3.3 Check the mapping rate to different Plant Genomes using bowtie

    $check_UBowtie.pl -i list -s SS -p CPU -d genome_bowtie_index -n genome_name

3.4 Parse report file for rRNA removing

    $parse_rRNA_rm_rpt.pl rRNA_rm.report > output

=================
4. Analysing Grafting RNASeq Datasets
=================

4.1 Aligning reads to scion sequence

    $tophat_pipeline_v2.pl -i list -d /path/scion_sequnece -s PE -l fr-unstranded -p 8 -m 1 -n 1    

4.2 Get unmapped reads
 
    $get_unmapped_read.pl -i list -s PE

4.3 Aligning unmapped reads to grated sequence

    $tophat_pipeline_v2.pl -i list_unmap -d /path/grated_sequnece -s PE -l fr-unstranded -p 8 -m 0 -n 0

4.4 Count the raw number for each gene & Get mapped reads and corresponding genes with fasta format

    $get_exp_raw.pl -i list_unmap -s PE -a grated_gene_position -r 

4.5 Filter the raw count (generated at step 4.4) without read support
    
    $filter_raw_count_graft.pl raw_count > output

