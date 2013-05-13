#!/usr/bin/perl

=head1

 Author: Yi Zheng
 Update: 02-26-2013

=cut

use strict;
use warnings;
use IO::File;
use Getopt::Long;

my $usage = qq'
Perl get_unmapped_read.pl -i list

 -i|list                (str)   input read list file (required)
 -o|output-prefix       (str)   prefix of out files (default = unmap)
 -s|sequencing-method   (str)   (required) (default = SS)
        PE (paired-end);
        SE (single-end);
        PS (paired strand-specific);
        SS (single strand-specific);
 -h|?|help                help info

Function:
Get unmapped reads after aligning RNA-Seq reads to genome using tophat.
The unmapped reads could be used for:
1. de novo assembly for novel genes.
2. check the contamination

';

my $help;
my ($read_list, $sequencing_method, $output_prefix);

GetOptions(
        "h|?|help"              => \$help,
        "i|list=s"              => \$read_list,
        "s|sequencing-method=s" => \$sequencing_method,
        "o|output-prefix=s"     => \$output_prefix
);

die $usage if $help;
die $usage unless $read_list;

$output_prefix ||= "unmap";
$sequencing_method ||= "SS";


# check sequencing method
if ($sequencing_method ne "SE" && $sequencing_method ne "SS" && $sequencing_method ne "PE" && $sequencing_method ne "PS" ) {
	die "Error at sequencing-method: $sequencing_method\n";
}
my $sequencing = $sequencing_method;

# get the file base on sequencing
my %list_read = check_list_file($read_list, $sequencing_method);

#################################################################
# main								#
#################################################################
foreach my $pre (sort keys %list_read)
{
	# get the reads and bam file from list hash
	my ($file, $file2, $plus_bam, $minus_bam, $bam);

        if ($sequencing_method eq "PE" || $sequencing_method eq "PS")
        {
                if ($sequencing_method eq "PE") {
                        ($file, $file2, $bam) = split(/\t/,$list_read{$pre});
                } else {
                        ($file, $file2, $plus_bam, $minus_bam) = split(/\t/,$list_read{$pre});
                }
		if ($file eq "NA" || $file2 eq "NA") { die "Error, clean read files $file and $file2 do not exist.\n"; }
        }
        else
        {
                if ($sequencing_method eq "SE") {
                        ($file, $bam) = split(/\t/, $list_read{$pre});
                } else {
                        ($file, $plus_bam, $minus_bam) = split(/\t/,$list_read{$pre});
                }
		if ($file eq "NA") { die "Error, clean read file $file does not exist.\n"; }
        }

	# mapped id to hash
	my %mapped;
        if ($sequencing_method eq "PS" || $sequencing_method eq "SS")
        {
                my $plus_sam = $pre."_plus.sam";
                my $minus_sam = $pre."_minus.sam";
                system("samtools view -h -o $plus_sam $plus_bam") && die "Error at samtools view -h -o $plus_sam $plus_bam\n";
                system("samtools view -h -o $minus_sam $minus_bam") && die "Error at samtools view -h -o $minus_sam $minus_bam\n";

                my $pfh = IO::File->new($plus_sam) || die "Can not open plus SAM file $plus_sam $!\n";
                while(my $line = <$pfh>)
                {
                        my @a = split(/\t/, $line);
                        unless ($line =~ m/^@/)
                        {
                                $mapped{$a[0]} = 1;
                        }
                }
                $pfh->close;

                my $mfh = IO::File->new($minus_sam) || die "Can not open minus file $minus_sam $!\n";
                while(my $line = <$mfh>)
                {
                        my @a = split(/\t/, $line);
                        unless ($line =~ m/^@/)
                        {
                                $mapped{$a[0]} = 1;
                        }
                }
                $mfh->close;

                unlink($plus_sam);
                unlink($minus_sam);
        }
	else
        {
                my $sam = $pre."_all.sam";
                system("samtools view -h -o $sam $bam") && die "Error at samtools view -h -o $sam $bam\n";

                my $sfh = IO::File->new($sam) || die "Can not open sam file $sam $!\n";
                while(my $line = <$sfh>)
                {
                        my @a = split(/\t/, $line);
                        unless ($line =~ m/^@/)
                        {
                                $mapped{$a[0]} = 1;
                        }
                }
                $sfh->close;

                #unlink($sam);
        }	

	# output unmapped reads
	my ($out_file1, $out_file2, $fh, $fh2, $format, $id, $id2);

	if ($sequencing_method eq "PS" || $sequencing_method eq "PE")
	{
		$out_file1  = $pre."_".$output_prefix."_1.fa";
		$out_file2 = $pre."_".$output_prefix."_2.fa";


		my $out1 = IO::File->new(">".$out_file1) || die "Can not open unmapped read file1 $!\n";
		my $out2 = IO::File->new(">".$out_file2) || die "Can not open unmapped read file2 $!\n";

		# open clean read file
		if ($file =~ m/\.gz/)
                {
                        open($fh, "<:gzip", $file) || die "Can not open read file $file $!\n";
                }
                else
                {
                        open($fh, $file) || die "Can not open read file $file $!\n";
                }
		if ($file2 =~ m/\.gz/)
                {
                        open($fh2, "<:gzip", $file2) || die "Can not open read file $file2 $!\n";
                }
                else
                {
                        open($fh2, $file2) || die "Can not open read file $file2 $!\n";
                }
	
		while($id = <$fh>)
                {
                        chomp($id);
                        if ($id =~ m/^>/) { $format = "fa";   $id =~ s/^>//; }
                        elsif ($id =~ m/^@/) { $format = "fq";$id =~ s/^@//; }
                        else { die "Error at read format $id\n"; }
                        my $seq = <$fh>; chomp($seq);
                        if ($format eq "fq") { <$fh>; <$fh>; }

			my @a = split(/\s+/, $id);
			#print $a[0]."\n"; die;
			unless( defined $mapped{$a[0]} ) { print $out1 ">".$id."\n".$seq."\n";  }
                }
                close($fh);

		$format = "";
		while($id2 = <$fh2>)
                {
                        chomp($id2);
                        if ($id2 =~ m/^>/) { $format = "fa";    $id2 =~ s/^>//; }
                        elsif ($id2 =~ m/^@/) { $format = "fq"; $id2 =~ s/^@//; }
                        else { die "Error at read format $id2\n"; }
                        my $seq2 = <$fh2>; chomp($seq2);
                        if ($format eq "fq") { <$fh2>; <$fh2>; }

			my @a = split(/\s+/, $id2);
                        unless( defined $mapped{$a[0]} ) { print $out2 ">".$id2."\n".$seq2."\n"; }
                }
                close($fh);

                $out1->close;
		$out2->close;
	}
	else
	{
		$out_file1 = $pre."_".$output_prefix.".fa";
		my $out = IO::File->new(">".$out_file1) || die "Can not open unmapped read file1 $!\n";

		# open clean read file
		if ($file =~ m/\.gz/)
                {
                        open($fh, "<:gzip", $file) || die "Can not open read file $file $!\n";
                }
                else
                {
                        open($fh, $file) || die "Can not open read file $file $!\n";
                }

		# read clean read file and output unmapped reads
		while($id = <$fh>)
                {
                        chomp($id);
			$format = "";
                        if ($id =~ m/^>/)    { $format = "fa"; $id =~ s/^>//; }
                        elsif ($id =~ m/^@/) { $format = "fq"; $id =~ s/^@//; }
                        else { die "Error at read format $id\n"; }

                        my $seq = <$fh>; chomp($seq);
                        if ($format eq "fq") { <$fh>; <$fh>; }
				
			my @a = split(/\s+/, $id);
			unless($mapped{$a[0]}) { print $out ">".$id."\n".$seq."\n"; }
                }
                close($fh);
		$out->close;
	}
}

#################################################################
# kentnf: subroutine                                            #
#################################################################

=head1 Sub:check_list_file

=cut
sub check_list_file
{
        my ($list_file, $sequencing_method) = @_;
        my %list_info;
        my @file_surfix = (".gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz", ".fq", ".fq.gz", ".fastq", ".fastq.gz");

        my ($reads, $bams);
        my $ls = IO::File->new($list_file) || die "Can not open list file $list_file $!\n";
        while(my $list = <$ls>)
        {
                chomp($list);

		# check if the read files are exist
		if ($sequencing_method eq "PS" || $sequencing_method eq "PE")
                {
                        my $file1 = $list."_1";
                        my $file2 = $list."_2";
                        if (-s $file1 && -s $file2)
                        {
                                $reads = $file1."\t".$file2;
                        }
                        else
                        {
                                foreach my $surfix (@file_surfix)
                                {
                                        if (-s $file1.$surfix && -s $file2.$surfix)
                                        {
                                                $reads = $file1.$surfix."\t".$file2.$surfix;
                                        }
                                }
                        }
                        unless ($reads) { $reads = "NA\tNA"; }
                }
                else
                {
                        if (-s $list)
                        {
                                $reads = $list;
                        }
                        else
                        {
                                foreach my $surfix (@file_surfix)
                                {
                                        if (-s $list.$surfix) {
                                                $reads = $list.$surfix;
                                        }
                                }
                        }
                        unless ($reads) { $reads = "NA"; }
                }

		# check bam file
		if ($sequencing_method eq "PS" || $sequencing_method eq "SS")
                {
                        my $plus_bam  = $list."_plus.bam";
                        my $minus_bam = $list."_minus.bam";
                        if (-s $plus_bam && -s $minus_bam) { $bams = $plus_bam."\t".$minus_bam; }
                }
                else
                {
                        my $bam = $list."_all.bam";
                        if (-s $bam) { $bams = $bam; }
                }

		# put reads and bam to hash;
		if($reads && $bams)
                {
                        print STDERR "Locate file for list: $list.\nRead File: $reads;\n BAM File: $bams.\n";
                        $list_info{$list} = $reads."\t".$bams;
                }
                else
                {
                        print STDERR "Error! Can not locate read and BAM for list: $list.\n";
                }
        }
        $ls->close;
        return %list_info;
}
