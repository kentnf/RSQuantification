#!/usr/bin/perl

use IO::File;
use PerlIO::gzip;

my $usage = qq'
Perl get_uniq_read.pl  read_file[fastq or fasta format, support gzip file]  Prefix_Read_ID

Function:
Parse the fasta or fastq file, then output the uniq fasta file with count

Example: 10930 is the number of seq in input file
>id000001-10930
ATGAAAAAAAAAAACCCCCCCCCCCCCCCCCC

';

my $read_file = shift || die $usage;
my $read_id_prefix = shift || 'ID';

my $out_file = $read_file;
$out_file =~ s/\..*$/_uniq\.fa/;

my %read = ();

my ($fh, $format, $total_read, $length);
if ($read_file =~ m/\.gz$/)
{
	$out_file =~ s/\.gz$//;
	open($fh, "<:gzip", $read_file) || die "Can not open read file $read_file $!\n";
}
else
{
	open($fh, $read_file) || die "Can not open read file $read_file $!\n";
}

while(my $id = <$fh>)
{
	chomp;
	$format = "";
	if ($id =~ m/^>/) { $format = 'fasta'; }
	elsif ($id =~ m/^@/) { $format = 'fastq'; }
	else { die "Can not identify sequence file format for ID $id\n";  }

	my $seq = <$fh>;
	chomp($seq);

	if ( defined $read{$seq} )
	{
		$read{$seq}++;
	}
	else
	{
		$read{$seq} = 1;
	}

	if ($format eq "fasta") { }
	else { <$fh>; <$fh>; }
	$total_read++;
}
$fh->close;

$length = length($total_read);

# sort the reads by num
my %sort_read_by_num;
foreach my $seq (sort keys %read)
{
	my $num = $read{$seq};
	if ($num > 1)
	{
		if (defined $sort_read_by_num{$num} )
        	{
                	$sort_read_by_num{$num}.="\t".$seq;
        	}
        	else
        	{
                	$sort_read_by_num{$num} = $seq;
        	}
		$read{$seq};
		delete $read{$seq};
	}
}

print scalar(keys(%sort_read_by_num))."\n";

# output the results
my $seq_num = 0;
open(OUT, ">".$out_file) || die "Can not open output file $out_file $!\n";
foreach my $num (sort { $b<=>$a } keys %sort_read_by_num)
{
	my @seq = split(/\t/, $sort_read_by_num{$num});
	foreach my $seq (@seq)
	{
		$seq_num++;
		my $seq_id = add_zero($seq_num, $length);
		$seq_id = $read_id_prefix.$seq_id;
		print OUT ">$seq_id-$num\n$seq\n";
	}
}

foreach my $seq (sort keys %read)
{
	$seq_num++;
	my $seq_id = add_zero($seq_num, $length);
	$seq_id = $read_id_prefix.$seq_id;
	print OUT ">$seq_id-1\n$seq\n";
}

close(OUT);

print "No. of read: $total_read\nNo. of uniq read: $seq_num\n";

#################################################################
# kentnf: subroutine						#
#################################################################

=head1

=cut
sub add_zero
{
	my ($id, $length) = @_;
	my $id_len = length($id);
	my $zero = "";
	for(my $i=0; $i<$length-$id_len; $i++)
	{
		$zero.="0";
	}
	my $return_id = $zero.$id;
	return $return_id;
}
