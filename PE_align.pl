#!/usr/bin/perl

=head

 make the PE reads alignments

=cut

use strict;
use warnings;
use IO::File;

my $usage = qq'\nusage: PE_align.pl read1 read2\n\n';

my $read1 = shift || die $usage;
my $read2 = shift || die $usage;

# read id to hash
my %gid1 = get_read_id($read1);
my %gid2 = get_read_seq($read2);
my %gid;
foreach my $id (sort keys %gid1) 
{  
	if (defined $gid2{$id})
	{
		$gid{$id} = 1;
	}
}

# PE alignment
my $output1 = $read1.".align";
my $output2 = $read2.".align";

my $out1 = IO::File->new(">".$output1) || die "Can not open output file1 $output1 $!\n";
my $out2 = IO::File->new(">".$output2) || die "Can not open output file1 $output2 $!\n";

my %ordered_reads; my $num = 0;

my $format;
my $fh1 = IO::File->new($read1) || die "Can not open read file1 $read1 $!\n";
while(<$fh1>)
{
	my $id1 = $_; chomp($id1);

	if	($id1 =~ m/^@/) { $format = 'fastq'; $id1 =~ s/^@//; }
	elsif 	($id1 =~ m/^>/) { $format = 'fasta'; $id1 =~ s/^>//; }
	else	{ die "Error in sequence format $id1\n"; }

	my @a = split(/\s+/, $id1, 2);
	my $seq = <$fh1>; chomp($seq);

	if (defined $gid{$a[0]})
	{
		$num++;
		$ordered_reads{$num} = $a[0];

		if ($format eq 'fastq')
		{
			my $id2 = <$fh1>; chomp($id2);
			my $qul = <$fh1>; chomp($qul);
			print $out1 "@"."$id1\n$seq\n$id2\n$qul\n";
		}
		else
		{
			print $out1 ">"."$id1\n$seq\n";
		}
	}
}
$fh1->close;

foreach my $n (sort {$a<=>$b} keys %ordered_reads )
{
	my $id = $ordered_reads{$n};
	print $out2 $gid2{$id};
}

$out1->close;
$out2->close;

#################################################################
# kentnf: subroutine						#
#################################################################
sub get_read_id
{
	my $file = shift;
	my %gid; my $format;
	my $fh = IO::File->new($file) || die "Can not open read file $file\n";
	while(<$fh>)
	{
		chomp;
		my $id = $_; 
		if 	($id =~ m/^@/ ) { $format = 'fastq'; $id =~ s/^@//; }
		elsif	($id =~ m/^>/ ) { $format = 'fasta'; $id =~ s/^>//; }
		else	{ die "Error in sequence format $id\n"; }

		my @a = split(/\s+/, $id);
		$gid{$a[0]} = 1;
		my $seq = <$fh>;

		if ($format eq 'fastq') { <$fh>; <$fh>; }
	}
	$fh->close;
	return %gid;
}

sub get_read_seq
{
        my $file = shift;
        my %gid; my $format;
        my $fh = IO::File->new($file) || die "Can not open read file $file\n";
        while(<$fh>)
        {
                chomp;
                my $id = $_;
                if      ($id =~ m/^@/ ) { $format = 'fastq'; $id =~ s/^@//; }
                elsif   ($id =~ m/^>/ ) { $format = 'fasta'; $id =~ s/^>//; }
                else    { die "Error in sequence format $id\n"; }

                my @a = split(/\s+/, $id);
                my $seq = <$fh>;

		my $reads;

                if ($format eq 'fastq') 
		{
			my $id2 = <$fh>;
			my $qul = <$fh>;
			$reads = "@".$id."\n".$seq.$id2.$qul;
		}
		else
		{
			$reads = ">".$id."\n".$seq;
		}

		$gid{$a[0]} = $reads;
        }
        $fh->close;
        return %gid;
}

