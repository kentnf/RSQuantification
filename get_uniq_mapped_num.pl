#!/usr/bin/perl

#===============================================================================
=head1 Name

 get_uniq_mapped_num.pl

=head1 Version 

 Author:   Yi Zheng
 Version:  1.0

=head1 Update

 2012-02-29
 1. support read gz file
 2. add total read num
 3. add mapped percentage
 4. add parameters for strand-specific and none strand-specific

 2012-03-09
 1. identify read format automatically.
<<<<<<< HEAD

 2013-02-23
 1. read file is not required, just bam is OK
=======
>>>>>>> 1f8cef8346c6584816b471a27f17b383d50a3df6
 
=head1 Description

	-i|list		 	(str)  input read list file (required)
	-o|output		(str)  out file name   (default = uniq_mapped_num)
	-s|sequencing-method	(str)  (required) (default = SS)
			PE (paired-end); 
			SE (single-end); 
			PS (paired strand-specific); 
			SS (single strand-specific);
	-h|?|help	 	 help info

=head1 Example

	For single end strand-specific sequencing
        perl get_uniq_mapped_num.pl -i read_list -s SS

=cut
#===============================================================================

use strict;
use warnings;
use IO::File;
use PerlIO::gzip;
use Getopt::Long;

my $help;
my ($list_file, $output, $sequencing_method);

GetOptions(
        "h|?|help"		=> \$help,
        "i|list=s"		=> \$list_file,
	"o|output=s"		=> \$output,
	"s|sequencing-method=s"	=> \$sequencing_method
);

die `pod2text $0` if $help;
die `pod2text $0` unless $list_file;

$output ||= "uniq_mapped_num";
$sequencing_method ||= "SS";

# check sequencing method
if ($sequencing_method ne "SE" && $sequencing_method ne "SS" && $sequencing_method ne "PE" && $sequencing_method ne "PS" ) {
        die "Error at sequencing-method: $sequencing_method\n";
}

# check list file
my %list_read = check_list_file($list_file, $sequencing_method);

my $out = IO::File->new(">".$output) || die "Can not open file $output $!\n";
<<<<<<< HEAD
print $out "#Sample\tclean\tmapped reads\t%mapped\n";
=======
>>>>>>> 1f8cef8346c6584816b471a27f17b383d50a3df6

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
        }
        else
        {
		if ($sequencing_method eq "SE") {
                	($file, $bam) = split(/\t/, $list_read{$pre});
		} else {
			($file, $plus_bam, $minus_bam) = split(/\t/,$list_read{$pre});
		}
        }
	
	# count the total number of read
	my ($fh, $id, $format, $total_read);
<<<<<<< HEAD
	if ($file eq "NA")
	{
		$total_read = "NA";
	}
	else
	{
		if ($file =~ m/\.gz/)
		{
			open($fh, "<:gzip", $file) || die "Can not open read file $file $!\n";
		}
		else
		{
			open($fh, $file) || die "Can not open read file $file $!\n";
		}

		while($id = <$fh>)
		{
			chomp($id);
			unless ($format)
			{
				if ($id =~ m/^>/) { $format = "fa"; }
				elsif ($id =~ m/^@/) { $format = "fq"; }
				else { die "Error at read format $id\n"; }
			}
			<$fh>;
			if ($format eq "fq") { <$fh>; <$fh>; }
			$total_read++;
		}
		close($fh);
	}
=======
	if ($file =~ m/\.gz/)
	{
		open($fh, "<:gzip", $file) || die "Can not open read file $file $!\n";
	}
	else
	{
		open($fh, $file) || die "Can not open read file $file $!\n";
	}

	while($id = <$fh>)
	{
		chomp($id);
		unless ($format)
		{
			if ($id =~ m/^>/) { $format = "fa"; }
			elsif ($id =~ m/^@/) { $format = "fq"; }
			else { die "Error at read format $id\n"; }
		}
		<$fh>;
		if ($format eq "fq") { <$fh>; <$fh>; }
		$total_read++;
	}
	close($fh);
>>>>>>> 1f8cef8346c6584816b471a27f17b383d50a3df6

	# step2. count the number of mapped read
	my ($uniq_mapped_read, $percentage);
		
	my %mapped;
	if ($sequencing_method eq "PS" || $sequencing_method eq "SS")	
	{
<<<<<<< HEAD
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
=======

		my $pfh; 
		open ($pfh, "samtools view $plus_bam | ") or die "Can not open plus BAM file $plus_bam with command \"samtools view $plus_bam\".\n$!\n"; 
		while(my $line = <$pfh>)
		{
			my @a = split(/\t/, $line);
			$mapped{$a[0]} = 1;
		}
		$pfh->close;

		my $mfh; 
		open ($mfh, "samtools view $minus_bam | ") or die "Can not open plus BAM file $minus_bam with \"samtools view $minus_bam\".\n$!\n"; 
		while(my $line = <$mfh>)
		{
			my @a = split(/\t/, $line);
                       	$mapped{$a[0]} = 1;
		}
		$mfh->close;
	
        }
	else
	{
		my $sfh; 
		open($sfh, "samtools view $bam") or die "Can not open bam file with command \"samtools view $bam\"\n$!\n"; 
>>>>>>> 1f8cef8346c6584816b471a27f17b383d50a3df6
                while(my $line = <$sfh>)
                {
                        my @a = split(/\t/, $line);
                        unless ($line =~ m/^@/)
                        {
                                $mapped{$a[0]} = 1;
                        }
                }
                $sfh->close;
<<<<<<< HEAD

                unlink($sam);
	}

        $uniq_mapped_read = scalar(keys(%mapped));

	if ($total_read eq "NA") { $percentage = "NA"; }
	else { $percentage =  sprintf("%.2f", ($uniq_mapped_read/$total_read)*100); }
        print $out "$pre\t$total_read\t$uniq_mapped_read\t$percentage\n";
=======
	}

        $uniq_mapped_read = scalar(keys(%mapped));
	$percentage =  sprintf("%.2f%%", ($uniq_mapped_read/$total_read)*100);
        print $out "$pre\t$uniq_mapped_read\t$total_read\t$percentage\n";
>>>>>>> 1f8cef8346c6584816b471a27f17b383d50a3df6
}
$out->close;


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
<<<<<<< HEAD
			unless ($reads) { $reads = "NA\tNA"; }
=======
>>>>>>> 1f8cef8346c6584816b471a27f17b383d50a3df6
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
<<<<<<< HEAD
			unless ($reads) { $reads = "NA"; }
=======
>>>>>>> 1f8cef8346c6584816b471a27f17b383d50a3df6
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
<<<<<<< HEAD
                        print STDERR "Error! Can not locate read and BAM for list: $list.\n";
=======
                        print STDERR "Error! Can not locate read or BAM for list: $list.\n";
>>>>>>> 1f8cef8346c6584816b471a27f17b383d50a3df6
                }
        }
	$ls->close;
        return %list_info;
}
