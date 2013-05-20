#!/usr/bin/perl

=head1 Name

 gff2gene_position.pl -- conert gff annotation to gene position file.

=cut

use strict;
use warnings;
use IO::File;
use Getopt::Long;

my $usage = qq'
usage: perl $0 input [options]

	h	help
	i=s  	input_file
	t=s	type
	n=s	name
	s=s	split

';

my ($help, $input_file, $type, $name, $split);

GetOptions(
	"h"	=> \$help,
	"i=s"	=> \$input_file,
	"t=s"	=> \$type,
	"n=s"	=> \$name,
	"s=s"	=> \$split
);

die $usage if $help;
die $usage unless $input_file;
$type ||= "mRNA";
$name ||= "ID";
$split ||= ";";

# main
my $fh = IO::File->new($input_file) || die "Can not open GTF/GFF file $input_file $!\n";
while(<$fh>)
{
	chomp;
	unless($_ =~ m/^#/)
	{
		my @a = split(/\t/, $_);
		if ($a[2] eq $type)
		{
			my ($chr, $start, $end, $strand, $annotation) = ($a[0], $a[3], $a[4], $a[6], $a[8]);
			my $gene_id = get_id($annotation, $type, $name, $split);
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
	my ($annotation, $name, $split) = @_;
	my @a = split(/\Q$split\E/, $annotation);
	my $gid;
	foreach my $a (@a)
	{
		my @b = split(/=/, $a);
		if ($b[0] eq $name)
		{
			$gid = $b[1];
		}
	}

	$gid =~ s/^\s+//ig;
	chomp($gid);
	return $gid;
}
