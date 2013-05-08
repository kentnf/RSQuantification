#!/usr/bin/perl

use strict;
use warnings;
use Statistics::Basic qw(:all nofill);

my $usage = "Perl this input > output\n";

my $file = shift || die $usage;

open(FH , $file) || die "Can not open input file: $file\n";

my $line_i = 0;

while(<FH>)
{
	chomp;

	@a = split(/\t/, $_);

	#########################################################
	# parse first line					#
	#########################################################
	$line_i++;

	if ($line_i == 1)
	{
		for(my $j=1; $j<@a; $j++)
		{
			my $ad_j = addzero($j, 2);

			$title{$ad_j} = $a[$j];
		}
	}

	#########################################################
	# parse data line					#
	#########################################################

	else
	{
		for(my $i=1; $i<@a; $i++)
		{
			my $ad_i = addzero($i, 2);

			if (defined $order{$ad_i})
			{
				$order{$ad_i}.="\t".$a[$i];
			}
			else
			{
				$order{$ad_i} = $a[$i];
			}
		}
	}
}


#################################################################
# produce correlation for every line data			#
#################################################################

#foreach my $a (sort keys %order) {  print $a."###".$order{$a}."\n";}


foreach my $ti (sort keys %title)
{
	print "\t".$title{$ti};
}
print "\n";

foreach my $oi (sort keys %order)
{
	print $title{$oi};

	my @array_i = split(/\t/, $order{$oi});

	foreach my $oj (sort keys %order)
	{
		my @array_j = split(/\t/, $order{$oj});

		my $correlation = corr([@array_i], [@array_j]);

		print "\t".$correlation;
	}

	print "\n";
}


sub addzero
{
	my ($num, $len) = @_;

	my $zero = "";

	for(my $i=0; $i<($len-length($num)); $i++)
	{
		$zero.="0";
	}

	return $zero.$num;
}

#my $correlation = corr( [1.2,2.3,4], [1.1,2.2,3.1] );
#print $correlation."\n";
