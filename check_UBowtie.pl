#!/usr/bin/perl

=head

 Check the mapping rate to different Plant Genomes using bowtie
 Author: Yi Zheng
 update: 03/10/2013

=cut

use strict;
use warnings;
use IO::File;
use Getopt::Long;

my $usage = q'
USAGE:
Perl check_UBowtie.pl -i list_fasta -s SequencingMethod (default = SS) -p CPU -d genome_bowtie_index -n genome_name

  SequencingMethod
  PE (paired-end);
  SE (single-end);
  PS (paired strand-specific);
  SS (single strand-specific);

* using this program to test how many read were mapped to different plant genome using bowtie

';

my $help;
my ($read_list, $sequencing_method, $cpu, $database_index, $database_name);

GetOptions(
        "h|?|help"              => \$help,
        "i|list=s"              => \$read_list,
        "s|sequencing-method=s" => \$sequencing_method,
	"d|database-index=s"	=> \$database_index,
	"n|database-name=s"	=> \$database_name,
	"p|cpu=i"		=> \$cpu
);

die $usage if $help;
die $usage unless $read_list;
die $usage unless $database_index;
die $usage unless $database_name;
$sequencing_method ||= "SS";
$cpu ||= 24;

#################################################################
# set plant genome sequence below				#
#################################################################
my %plantgenome = (
#'cucumber' => '/home/fei/kentnf/pub/PlantGenome/cucumber/cucumber_genome_v2',
#'tomato' => '/home/fei/kentnf/pub/PlantGenome/tomato/tomato_genome_240',
#'watermelon' => '/home/fei/kentnf/pub/PlantGenome/watermelon/watermelon_genome',
#'apple' => '/home/fei/kentnf/pub/PlantGenomes/apple/apple_genome',
#'grape' => '/home/fei/kentnf/pub/PlantGenomes/grape/grape_genome',
#'arabidopsis' => '/home/fei/kentnf/pub/PlantGenomes/arabidopsis/Arabidopsis_genome',
#'NBeth' => '/home/fei/kentnf/pub/PlantGenomes/NBeth/NBeth_genome',
#'potato' => '/home/fei/kentnf/pub/PlantGenomes/potato/potato_genome'
);

# check sequencing method
if ($sequencing_method ne "SE" && $sequencing_method ne "SS" && $sequencing_method ne "PE" && $sequencing_method ne "PS" ) {
        die "Error at sequencing-method: $sequencing_method\n";
}

# check list file
my %list_read = check_list_file($read_list, $sequencing_method);
my %return_list;

# check the genome index
my @index_files = (
        $database_index.".1.ebwt",
        $database_index.".2.ebwt",
        $database_index.".3.ebwt",
        $database_index.".4.ebwt",
        $database_index.".rev.1.ebwt",
        $database_index.".rev.2.ebwt"
);
foreach my $index_file (@index_files) {
        unless(-s $index_file) { die "Error! $index_file not exist or have no info\n"; }
}
$plantgenome{$database_name} = $database_index;

# main 
foreach my $plant (sort keys %plantgenome )
{
	my $plant_genome_index = $plantgenome{$plant};

	foreach my $list (sort keys %list_read)
	{
		if ($sequencing_method eq "SS" || $sequencing_method eq "SE")
		{
			my $read_file = $list_read{$list};
			my $mapped_file = $list."_".$plant."_mapped.fa";
			my $unmap_file = $list."_".$plant."_unmap.fa";
			my $cmd = "bowtie -v 0 -k 1 -p $cpu --un $unmap_file -f $plant_genome_index $read_file $mapped_file";
			print $cmd."\n";
			system($cmd) && die "Error at $cmd\n";
			$return_list{$list} = $unmap_file;
		}
		else
		{
			die;
			#my @read = split(/\t/, $list_read{$list});			
			#$cmd = "bowtie -v 0 -k 1 -p $cpu --un $list -f $plant_genome -1 $read[0] -2 $read[1] MMMM";
			#print $cmd."\n";
			#system($cmd) && die "Error at $cmd\n";
			#$return_list{$list} = $unmap_file1."\t".$unmap_file2;
		}
	}

	# generate new list using unmapped reads
	%list_read = %return_list;
	%return_list = ();
}

#################################################################
# kentnf: subroutine			 			#
#################################################################

=head1 Sub:check_list_file

=cut
sub check_list_file
{
        my ($list_file, $sequencing_method) = @_;
        my %list_info;
	my @file_surfix = (".fasta", ".fa");
        #my @file_surfix = (".gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz", ".fq", ".fq.gz", ".fastq", ".fastq.gz");

        my ($reads, $bams);
        my $ls = IO::File->new($list_file) || die "Can not open list file $list_file $!\n";
        while(my $list = <$ls>)
        {
                chomp($list);
                # check if the read files are exist
                if ($sequencing_method eq "PS" || $sequencing_method eq "PE")
                {
                        my $file1 = $list."_1";
                        my $file2 = $list."_2";
                        if (-s $file1 && -s $file2)
                        {
                                $reads = $file1."\t".$file2;
                        }
                        else
                        {
                                foreach my $surfix (@file_surfix)
                                {
                                        if (-s $file1.$surfix && -s $file2.$surfix)
                                        {
                                                $reads = $file1.$surfix."\t".$file2.$surfix;
                                        }
                                }
                        }
                }
                else
                {
                        if (-s $list)
                        {
                                $reads = $list;
                        }
                        else
                        {
                                foreach my $surfix (@file_surfix)
                                {
                                        if (-s $list.$surfix) {
                                                $reads = $list.$surfix;
                                        }
                                }
                        }
                }

                # put reads and bam to hash;
                if($reads )
                {
                        print STDERR "Locate file for list: $list.\nRead File: $reads;\n";
                        $list_info{$list} = $reads;
                }
                else
                {
                        print STDERR "Error! Can not locate read or BAM for list: $list.\n";
                }
        }
        $ls->close;
        return %list_info;
}
