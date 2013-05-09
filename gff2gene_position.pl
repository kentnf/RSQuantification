#!/usr/bin/perl

=head1 Name

 gff2gene_position.pl -- conert gff annotation to gene position file.

=cut

use strict;
use warnings;
use IO::File;

my $usage = qq'
usage: perl $0 input > output

* input is gff file, the gff file must have gene label
';

my $file = shift || die $usage;

my $fh = IO::File->new($file) || die "Can not open GTF file $file $!\n";
while(<$fh>)
{
	chomp;
	unless($_ =~ m/^#/)
	{
		my @a = split(/\t/, $_);
		if ($a[2] eq 'gene')
		{
			my ($chr, $start, $end, $strand) = ($a[0], $a[3], $a[4], $a[6]);
			my $gene_id = get_gene_id($a[8]);
			print "$chr\t$start\t$end\t$gene_id\t$strand\n";
		}
	}
}
$fh->close;

#################################################################
# kentnf: subroutine						#
#################################################################
=head2 get_gene_id

=cut
sub get_gene_id
{
	my $annotation = shift;
	my @a = split(/;/, $annotation);
	my $gid;
	foreach my $a (@a)
	{
		if ($a =~ m/^ID=/)
		{
			$gid = $a;
		}
	}
	$gid =~ s/\s+//ig;
	$gid =~ s/ID=//ig;
	return $gid;
}
