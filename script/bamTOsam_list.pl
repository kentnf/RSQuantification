#!/usr/bin/perl

use IO::File;

my $usage = q'
Perl bamTOsam.pl  list_of_bam

* convert bam file to sam file.
* this script relies on samtools, please install samtools first.
';

my $list = shift || die $usage;

my $fh = IO::File->new($list) || die "Can not open file $list\n";

while(<$fh>)
{
	chomp;
	$bam = $_;
	$sam = $bam;
	$sam =~ s/\.bam/\.sam/;

	system("samtools view -h -o $sam $bam") && die "Error at samtools view -h -o $sam $bam\n";

	print "$bam is convert to $sam successful\n";
}
$fh->close;
