#!/usr/bin/perl

#===============================================================~
=head1 Name

  get_exp_rpkm.pl

=head1 Version

 Author:    Yi Zheng
 version:   1.0

=head1 Update

 2013-12-18
 1. add option for libraries (cleaned reads, mapped reads)

 2012-03-05
 1. add help info and Getopt::Long

 2013-02-22
 1. using raw count and uniq_read_num to get RPKM

=head1 Descripton

 -x|gene-length	(str)	gene length
 -e|exp-raw	(str)	raw count
 -u|uniq-read	(str)	uniq mapped reads num
 -o|output	(str)	output (default = exp_rpkm)
 -b|library	(str)	for library selection : map/clean (default = map) 

=head1 Example
  
 Perl get_exp_rpkm.pl -e exp_sense_raw -u uniq_mapped_read_num -x tomato_gene_length

=head1 Input File Format

 gene length:
 gene_id \t length

 adjusted expression:
 gene \t sample1 \t sample2 ... sampleN
 libsize \t size1 \t size2 ... sizeN
 gene1 \t 0 \t 2 ... 10
 gene2 \t 5 \t 8 ... 40
 ..... 
 geneN \t 0 \t 6 ... 30

=cut
#================================================================

use strict;
use warnings;
use IO::File;
use Getopt::Long;

my $help;
my ($gene_length, $exp_raw, $uniq_read, $output);

my $library_type = 'map';

GetOptions(
	"h|?|help"		=> \$help,
	"x|gene-length=s"	=> \$gene_length,
	"u|uniq-read=s"		=> \$uniq_read,
	"e|exp-raw=s"		=> \$exp_raw,
	"b|library=s"		=> \$library_type,
	"o|output=s"		=> \$output
);

die `pod2text $0` if $help;
die `pod2text $0` unless $gene_length;
die `pod2text $0` unless $exp_raw;
die `pod2text $0` unless $uniq_read;
if ($library_type eq "map" || $library_type eq "clean") { } else { die die `pod2text $0`; }

$output ||= "exp_rpkm";
#================================================================
# put gene length to hash					#
#================================================================
my %gene_length;

my $afh = IO::File->new($gene_length) || die "Can not open gene length file $gene_length $!\n";
while(<$afh>)
{
	chomp;
	my @a = split(/\t/, $_);
	$gene_length{$a[0]} = $a[1];
}
$afh->close;

#================================================================
# Get libsize base on uniq mapped reads				#
#================================================================
my %lib_size;
my $ls = IO::File->new($uniq_read) || die "Can not open file $uniq_read $!\n";
while(<$ls>)
{
	chomp;
	unless($_ =~ m/^#/)
	{
		my @a = split(/\t/, $_);
		if ($library_type eq 'map') { $lib_size{$a[0]} = $a[2]; }
		else { $lib_size{$a[0]} = $a[1]; }
	}
}
$ls->close;

#================================================================
# compute RPKM using adjust raw expression, gene length, and library size
#================================================================
my ($length, $lib_size, $rpkm);

my $out = IO::File->new(">".$output) || die "Can not open output file $output $!\n";
my $efh = IO::File->new($exp_raw) || die "Can not open exp_adjust file $exp_raw $!\n";
my $title = <$efh>;
print $out $title;

chomp($title);
my @title = split(/\t/, $title);

while(<$efh>)
{
	chomp;
	my @a = split(/\t/, $_);
		
	my $this = "";
	my $print = 0;

	if (defined $gene_length{$a[0]})
	{
		$length = $gene_length{$a[0]};
		$this.= $a[0];
	}
	else
	{
		die "Error, no gene length for geneID: $a[0]\n$_\n";
	}

	for(my $j=1; $j<@a; $j++)
	{
		$lib_size = $lib_size{$title[$j]};
		die "Error at libsize : $lib_size for title: $title[$j]\n" unless $lib_size > 0;
		$rpkm = ($a[$j] * 1000 * 1000000) / ($length * $lib_size);
		$rpkm = sprintf("%.2f", $rpkm);
		$this.= "\t".$rpkm;
		if ($rpkm > 0) {$print = 1;}
	}

	print $out $this."\n";
}
$afh->close;
$out->close;

#================================================================
# kentnf : subroutine
#================================================================
