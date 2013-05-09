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

=head1 Descripton

 -x|gene-length	(str)	gene length
 -e|exp-adjust	(str)	adjusted expression
 -o|output	(str)	output (default = exp_rpkm)

=head1 Example
  
 Perl get_exp_rpkm.pl -e exp_sense_adjust -x tomato_gene_length

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
my ($method, $gene_length, $exp_adjust, $output);

GetOptions(
	"h|?|help"		=> \$help,
	"m|method=s"		=> \$method,
	"x|gene-length=s"	=> \$gene_length,
	"e|exp-adjust=s"	=> \$exp_adjust,
	"o|output=s"		=> \$output
);

die `pod2text $0` if $help;
die `pod2text $0` unless $gene_length;
die `pod2text $0` unless $exp_adjust; 
die `pod2text $0` unless $method;

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
# compute RPKM using adjust raw expression, gene length, and library size
#================================================================
my ($length, $lib_size, $rpkm, $rpm);

my $out = IO::File->new(">".$output) || die "Can not open output file $output $!\n";
my $efh = IO::File->new($exp_adjust) || die "Can not open exp_adjust file $exp_adjust $!\n";
my $title = <$efh>;
print $out $title;
my $lib_size_line = <$efh>;
chomp($lib_size_line);
my @lib_size = split(/\t/, $lib_size_line);
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
		$lib_size = $lib_size[$j];

		# for RPKM
		#$rpkm = ($a[$j] * 1000 * 1000000) / ($length * $lib_size);
		#$rpkm = sprintf("%.2f", $rpkm);
		#$this.= "\t".$rpkm;
		#if ($rpkm > 0) {$print = 1;}

		# for RPM
		$rpm = ($a[$j] * 1000000) / $lib_size;
		$rpm = sprintf("%.2f", $rpm);
		$this.= "\t".$rpm;
	}

	print $out $this."\n";
}
$afh->close;
$out->close;

#================================================================
# kentnf : subroutine
#================================================================
