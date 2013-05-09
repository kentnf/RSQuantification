#!/usr/bin/perl

#===============================================================~
=head1 Name

  get_exp_rpkm.pl

=head1 Version

 Author:    Yi Zheng
 version:   1.0

=head1 Update

 2012-03-05
 1. add help info and Getopt::Long

<<<<<<< HEAD
 2013-02-22
 1. using raw count and uniq_read_num to get RPKM

=head1 Descripton

 -x|gene-length	(str)	gene length
 -e|exp-raw	(str)	raw count
 -u|uniq-read	(str)	uniq mapped reads num
=======
=head1 Descripton

 -x|gene-length	(str)	gene length
 -e|exp-adjust	(str)	adjusted expression
>>>>>>> 1f8cef8346c6584816b471a27f17b383d50a3df6
 -o|output	(str)	output (default = exp_rpkm)

=head1 Example
  
<<<<<<< HEAD
 Perl get_exp_rpkm.pl -e exp_sense_raw -u uniq_mapped_read_num -x tomato_gene_length
=======
 Perl get_exp_rpkm.pl -e exp_sense_adjust -x tomato_gene_length
>>>>>>> 1f8cef8346c6584816b471a27f17b383d50a3df6

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
<<<<<<< HEAD
my ($gene_length, $exp_raw, $uniq_read, $output);
=======
my ($gene_length, $exp_adjust, $output);
>>>>>>> 1f8cef8346c6584816b471a27f17b383d50a3df6

GetOptions(
	"h|?|help"		=> \$help,
	"x|gene-length=s"	=> \$gene_length,
<<<<<<< HEAD
	"u|uniq-read=s"		=> \$uniq_read,
	"e|exp-raw=s"		=> \$exp_raw,
=======
	"e|exp-adjust=s"	=> \$exp_adjust,
>>>>>>> 1f8cef8346c6584816b471a27f17b383d50a3df6
	"o|output=s"		=> \$output
);

die `pod2text $0` if $help;
die `pod2text $0` unless $gene_length;
<<<<<<< HEAD
die `pod2text $0` unless $exp_raw;
die `pod2text $0` unless $uniq_read;
=======
die `pod2text $0` unless $exp_adjust; 
>>>>>>> 1f8cef8346c6584816b471a27f17b383d50a3df6

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
<<<<<<< HEAD
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
		$lib_size{$a[0]} = $a[2];
	}
}
$ls->close;

#================================================================
=======
>>>>>>> 1f8cef8346c6584816b471a27f17b383d50a3df6
# compute RPKM using adjust raw expression, gene length, and library size
#================================================================
my ($length, $lib_size, $rpkm);

my $out = IO::File->new(">".$output) || die "Can not open output file $output $!\n";
<<<<<<< HEAD
my $efh = IO::File->new($exp_raw) || die "Can not open exp_adjust file $exp_raw $!\n";
my $title = <$efh>;
print $out $title;

chomp($title);
my @title = split(/\t/, $title);

=======
my $efh = IO::File->new($exp_adjust) || die "Can not open exp_adjust file $exp_adjust $!\n";
my $title = <$efh>;
print $out $title;
my $lib_size_line = <$efh>;
chomp($lib_size_line);
my @lib_size = split(/\t/, $lib_size_line);
>>>>>>> 1f8cef8346c6584816b471a27f17b383d50a3df6
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
		die "Error, no gene length for $a[0]\n";
	}

	for(my $j=1; $j<@a; $j++)
	{
<<<<<<< HEAD
		#$lib_size = $lib_size[$j];
		$lib_size = $lib_size{$title[$j]};
=======
		$lib_size = $lib_size[$j];
>>>>>>> 1f8cef8346c6584816b471a27f17b383d50a3df6
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
