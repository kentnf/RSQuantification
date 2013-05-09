#!/usr/bin/perl

use IO::File;

my $file = shift || die;


my $fh = IO::File->new($file) || die "Can not open GTF file $file $!\n";
while(<$fh>)
{
	chomp;
	unless($_ =~ m/^#/)
	{
		my @a = split(/\t/, $_);
<<<<<<< HEAD
		if ($a[2] eq 'gene')
=======
		if ($a[2] eq 'mRNA')
>>>>>>> 1f8cef8346c6584816b471a27f17b383d50a3df6
		{
			my ($chr, $start, $end, $strand) = ($a[0], $a[3], $a[4], $a[6]);
			my $gene_id = get_gene_id($a[8]);
			print "$chr\t$start\t$end\t$gene_id\t$strand\n";
		}
	}
}
$fh->close;

# subroutine
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

