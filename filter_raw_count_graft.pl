#!/usr/bin/perl

use strict;
use warnings;
use IO::File;

my $usage = qq'\nusage: $0 raw_cout > output\n* this program is used for filtering raw count for graft dataset\n\n';
my $raw_count = shift || die $usage;

my $fh = IO::File->new($raw_count) || die "Can not open raw count file $raw_count\n";
my $title = <$fh>;
print $title;

my $total;
while(<$fh>)
{
	chomp;
	my @a = split(/\t/, $_);
	$total = 0;
	for(my $i=1; $i<@a; $i++) { $total = $total + $a[$i]; }
	if ( $total > 0 ) { print $_."\n"; }
}
$fh->close;
