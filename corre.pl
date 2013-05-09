#!/usr/bin/perl

=head1 Name

 corre.pl -- get correlation for each sample

 Author: Yi Zheng
 Email: yz357 (at) cornell.edu 

=cut

use strict;
use warnings;
use IO::File;
use Statistics::Basic qw(:all nofill);

my $usage = qq'
usage: perl $0 input > output

* input is expression of RPKM, output is correlation
';

my $file = shift || die $usage;

my %sample_exp;

my $fh = IO::File->new($file) || die "Can not open input file: $file\n";
# parse first title line
my $title = <$fh>; chomp($title);
my @title = split(/\t/, $title);

while(<$fh>)
{
	chomp;
	my @a = split(/\t/, $_);
	for(my $i=1; $i<@a; $i++)
	{
		my $tt = $title[$i];
		if (defined $sample_exp{$tt})	{ $sample_exp{$tt}.="\t".$a[$i]; }
		else				{ $sample_exp{$tt} = $a[$i]; }
	}
}
$fh->close;

#################################################################
# produce correlation for every line data			#
#################################################################

print "Sample";
shift @title;
foreach my $tt ( @title ) { print "\t".$tt; }
print "\n";

foreach my $tt_i ( @title )
{
	print $tt_i;
	my @exp_i = split(/\t/, $sample_exp{$tt_i});

	foreach my $tt_j ( @title )
	{
		my @exp_j = split(/\t/, $sample_exp{$tt_j});
		my $correlation = corr([@exp_i], [@exp_j]);
		print "\t".$correlation;
	}
	print "\n";
}

