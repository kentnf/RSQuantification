#!/usr/bin/perl

#================================================================
=head1 Name

 tophat_pipeline_v2.pl

=head1 Version

 Author:  Yi Zheng
 Version: 2.0

=head1 Update
 
 2012/02/25 
 1. add parameter segment length
 2. add PS (paired strand-specific) sequencing method
 please test these two new function ASAP 

 2012/03/07
 1. support gz file
 2. add -N parameter for v1.4.0
 3. use Getopt::Long modular

 2012/03/18
 1. identify fa and fq format automatically

 2012/04/17
 1. fix two bugs, a) get correct id from fa or fq file; b) parse alignemtn with mismatch

 2012/09/9
 1. fix bug for parameter cycle 
 2. change script to tophat_pipeline_v1.pl for tophat v1
 3. run cufflinks in bin first, than find cufflinks in path
 4. fix bug for the ID of fq sequence

=head1 Example
 
=cut
#use strict;
#use warnings;
use FindBin;
use IO::File;
use PerlIO::gzip;
use Bio::SeqIO;
use Getopt::Long;

BEGIN { $ENV{PATH}="${FindBin::RealBin}:$ENV{PATH}"; } 

my $usage = qq'
tophat_pipeline_v2.pl

 Author:  Yi Zheng
 Version: 2.0

USAGE: tophat_pipeline_v2.pl -i list -s SS -d database -l fr-firststrand [options]

 -i  [String]   Name of the input list of read files (required)
 -s  [String]   Sequencing method for read files: (required)
        PE (paired-end);
        SE (single-end);
        PS (paired strand-specific);
        SS (single strand-specific);
 -d  [String]   Indexed database of genome sequences (required)
 -l  [String]   Library type, fr-unstranded, fr-firststrand, fr-secondstrand (required)
 -p  [Integer]  number of CPUs used for megablast clustering (default = 1)
 -r  [Integer]  mate-inner-dist (default = 100)
 -v  [Integer]  mate-std-dev (default = 20)
 -m  [Integer]  segment mismatches (default = 0)
 -n  [Integer]  initial read mismatches (default = 0)
 -b  [Integer]  setment length (default : half of read length; 17~25)
 -c  [Integer]  cycle (default = 0, range: 0,1,2)

 -a     (str)  gene and position file (optional)
 -x     (str)  gene length file (optional)

';

my $help;
my ($list_file, $sequencing_method, $database_index, $lib_type, $seg_mismatch, $seg_length, $read_mismatch, 
$mate_inner_dist, $mate_std_dev, $cycle, $gene_position, $gene_length);

GetOptions(
	"h|?|help"			=> \$help,
	"i|list=s"			=> \$list_file,
	"s|sequencing-method=s"		=> \$sequencing_method,
	"l|library-type=s"		=> \$library_type,

	"m|segment-mismatch=i"		=> \$seg_mismatch,
	"b|segment-length=i"		=> \$seg_length,
	"n|read-mismatch=i"		=> \$read_mismatch,

	"d|database-index=s"		=> \$database_index,
	"a|gene-position=s"		=> \$gene_position,
	"x|gene-length=s"		=> \$gene_length,

	"r|mate-inner-distance=i"	=> \$mate_inner_dist,
	"v|mate-std-dev=i"		=> \$mate_std_dev,

	"p|cpu=i"			=> \$cpu,
	"c|cycle=i"			=> \$cycle
);

print STDERR "\nTophat pipeline begin at:\t".`date`."\n";
print STDERR "Checking input parameters ......\n";

# check parameters.
# 1. check required parameters
# 2. set default parameters
# 3. check and correct the parameters
die $usage if $help;
die $usage unless $list_file;
die $usage unless $database_index;
$sequencing_method ||= "SS";
$library_type ||= "fr-firststrand";
$seg_mismatch ||= "0";
$read_mismatch ||= "0";
$mate_inner_dist ||= "100";
$mate_std_dev ||= "20";
$cpu ||= "8";
$cycle ||= "0";

if ($library_type ne "fr-unstranded" && $library_type ne "fr-firststrand" && $library_type ne "fr-secondstrand") {
	die "Error at library-type: $library_type\t".$usage."\n";
}
if ($sequencing_method ne "SE" && $sequencing_method ne "SS" && $sequencing_method ne "PE" && $sequencing_method ne "PS" ) {
	die "Error at sequencing-method: $sequencing_method\t".$usage."\n";
}
if ($cpu < 1) { $cpu = 1;  print STDERR "Change the number of used CPU to 1 (mininum)\n"; }
if ($cpu > 20) {$cpu = 20; print STDERR "Change the number of used CPU to 20 (maxinum)\n"; }
# if ($seg_mismatch < $read_mismatch) { $seg_mismatch = $read_mismatch; print STDERR "Change segment mismatch to $read_mismatch base on read mismatch\n"; }
if ($seg_mismatch < 0) { $seg_mismatch = 0;  print STDERR "Change segment mismatch to 0 (mininum)\n"; }
if ($seg_mismatch > 5) { $seg_mismatch = 5;  print STDERR "Change segment mismatch to 5 (maxinum)\n"; }
if ($read_mismatch < 0) { $read_mismatch = 0;  print STDERR "Change read mismatch to 0 (mininum)\n"; }
if ($read_mismatch > 5) { $read_mismatch = 5;  print STDERR "Change read mismatch to 5 (maxinum)\n"; }
if ($cycle < 0) {$cycle = 0; print STDERR "Change run cycle to 0 (mininum)\n"; }
if ($cycle > 2) {$cycle = 2; print STDERR "Change run cycle to 2 (maxinum)\n"; }
if ($seg_length) {
	if ($seg_length < 17) { $seg_length = 17; print STDERR "Change segment length to 17 (mininum)\n";}
	if ($seg_length > 25) { $seg_length = 25; print STDERR "Change segment length to 25 (maxinum)\n";}
}
# check input files
my @index_files = (
	$database_index.".1.bt2",
	$database_index.".2.bt2",
	$database_index.".3.bt2",
	$database_index.".4.bt2",
	$database_index.".rev.1.bt2",
	$database_index.".rev.2.bt2"
);

foreach $index_file (@index_files) {
	unless(-s $index_file) { die "Error! $index_file not exist or have no info\n"; }
}

# check parameters mate-inner-distance and mate-std-dev 
if ($sequencing_method eq "PE" || $sequencing_method eq "PS") {
	# find the range of mate_inner_dist and mate_std_dev
	if ($mate_inner_dist && $mate_std_dev) {}
	else { die "Error! Please use suitable mate_inner_dist and mate_std_dev for paired datasets\n".$usage."\n";}
}

# check library-type of strand-specific datasets
my $split_plus_minus = 0;
if ($sequencing_method eq "PS" || $sequencing_method eq "SS") {
	$split_plus_minus = 1;
	if ($library_type ne "fr-firststrand" && $library_type ne "fr-secondstrand") {
		die "Error! Please use suitable library type for your strand-specific dataset\n".$usage;
	}
} else {
	if ($library_type ne "fr-unstranded") {
		$library_type = "fr-unstranded";
		print STDERR "library type is changed to fr-unstranded on base of your sequencing method (PE or SE)\n";
	}
}

# check input list files
# check paired file or single file base on sequencing method
# put corrected files and corresponding list to hash
my %list = check_list_file($list_file, $sequencing_method);

#################################################################
# set basic parameters						#
#################################################################
$library_type = "--library-type $library_type";
$cpu = "-p $cpu";
$mate_inner_dist = "-r $mate_inner_dist";
$mate_std_dev = "--mate-std-dev $mate_std_dev";
$seg_mismatch = "--segment-mismatches $seg_mismatch";
$read_mismatch = "--read-mismatches $read_mismatch";

# add other parameters to basic for tophat (optional)			#
my $other_parameters = "";
# example 
# other_parameters = "-G iTAGv2.1.gtf";

my $basic_parameters;

if ($sequencing_method eq "PE" || $sequencing_method eq "PS")
{
	$basic_parameters = "--no-convert-bam $library_type $mate_inner_dist $mate_std_dev $cpu $seg_mismatch $read_mismatch $other_parameters";
}
else
{
	$basic_parameters = "--no-convert-bam $library_type $cpu $seg_mismatch $read_mismatch $other_parameters";
}

print $basic_parameters."\n";


# check tophat bin
my $tophat_bin;
if (-s "${FindBin::RealBin}/bin/tophat-2.0.4.Linux_x86_64/tophat2" )
{
	$tophat_bin = "${FindBin::RealBin}/bin/tophat-2.0.4.Linux_x86_64/tophat2";
}
else
{
	$tophat_bin = "tophat2";
}

print "\ntopaht bin locate at:\n".$tophat_bin."\n";

#================================================================
# main								
#================================================================
foreach my $list (sort keys %list)
{
	# step 1: prepare init files and store read to hash
	my %seq_hash = (); my %seq_hash1 = (); my %seq_hash2 = ();
	my ($file, $file1, $file2, $head, $pre_seg_len, $pre_seg_len1, $pre_seg_len2, $seq_hash, $seq_hash1, $seq_hash2);
	$head = $list."_head";

	if ($sequencing_method eq "PE" || $sequencing_method eq "PS")
	{
		($file1,$file2) = split(/\t/,$list{$list});
		($seq_hash1, $pre_seg_len1) = seq_to_hash($file1, $cycle);
		($seq_hash2, $pre_seg_len2) = seq_to_hash($file2, $cycle);
		%seq_hash1 = %$seq_hash1; %seq_hash2 = %$seq_hash2;
		if ($pre_seg_len1 < $pre_seg_len2) { $pre_seg_len = $pre_seg_len1; }
		else {$pre_seg_len = $pre_seg_len2; }
	}
	else
	{
		$file = $list{$list};
		($seq_hash, $pre_seg_len) = seq_to_hash($file, $cycle);
		%seq_hash = %$seq_hash;
	}

	if ($seg_length) { $pre_seg_len = $seg_length;}

	# set segment length to 17-25
	if ($pre_seg_len < 17) { $pre_seg_len = 17; print STDERR "Change segment length to 17 (mininum)\n";}
	if ($pre_seg_len > 25) { $pre_seg_len = 25; print STDERR "Change segment length to 25 (maxinum)\n";}

	if ($cycle == 0) {
		print "Seq Hash 1: ".scalar(keys(%seq_hash1))."\n";
		print "Seq Hash 2: ".scalar(keys(%seq_hash2))."\n";
		print "Seq Hash: ".scalar(keys(%seq_hash))."\n";
		print "Segment length: ".$pre_seg_len."\n";
	}

	#########################################################
	# start cycles for tophat				#
	#########################################################
	for(my $i=0; $i<=$cycle; $i++)
	{
		# step 2.1: align file1 file2 using tophat
		my $output_folder = $list."_m".$i;
		my $cmd_align;
		if ($sequencing_method eq "PE" || $sequencing_method eq "PS")
		{
			$cmd_align = "$tophat_bin -o $output_folder $basic_parameters --segment-length $pre_seg_len $database_index $file1 $file2";
		}
		else
		{
			$cmd_align = "$tophat_bin -o $output_folder $basic_parameters --segment-length $pre_seg_len $database_index $file";
		}
		print $cmd_align."\n"; 
		#die;
		system($cmd_align) && die "Error in commmand $cmd_align\n";
	
		# step 2.2 convert bam file to sam file
		my $output_sam = $output_folder."/accepted_hits.sam";

		# step 2.3 filter sam file mismatch (paired-end filter)
		my $output_sam_m = $output_folder."/mismatch.sam";
		my $temp_sam_m = $output_folder."/temp.sam";		
		my %select_read = ();
	       
		if ($sequencing_method eq "PE" || $sequencing_method eq "PS")
		{
			%select_read = paired_end_filter($output_sam, $output_sam_m, $temp_sam_m, $read_mismatch, $cycle);
		}
		else
		{
			%select_read = single_end_filter($output_sam, $output_sam_m, $read_mismatch, $cycle);
		}

		if ($cycle == 0)
		{
			print "Select Read: ".scalar(keys(%select_read))."\n";
		}
		
		# step 2.4 get unmaped read and save in new infile
		if ($i < $cycle)
		{
			my $j = $i+1;

			if ($sequencing_method eq "PE" || $sequencing_method eq "PS")
			{
				my ($num_unmapped1, $num_unmapped2);

				my $tmp_seq1 = $list."_read".$j."_1.gz";
				my $tmp_seq2 = $list."_read".$j."_2.gz";

				foreach my $select_read_id (sort keys %select_read)
				{
					delete $seq_hash1{$select_read_id};
					delete $seq_hash2{$select_read_id};
				}

				open(TFH1, ">:gzip", $tmp_seq1) || die "Can not open temp seq1 (unmapped) file: $tmp_seq1\n";
				foreach my $id (sort keys %seq_hash1)
				{
					print TFH1 ">".$id."\n".$seq_hash1{$id}."\n";
					$num_unmapped1++;
				}
				close(TFH1);

				open(TFH2, ">:gzip", $tmp_seq2) || die "Can not open temp seq2 (unmapped) file: $tmp_seq2\n";
				foreach my $id (sort keys %seq_hash2)
				{
					print TFH2 ">".$id."\n".$seq_hash2{$id}."\n";
					$num_unmapped2++;
				}
				close(TFH2);

				print "No. of unmapped reads in seq1: $num_unmapped1\n";
				print "No. of unmapped reads in seq2: $num_unmapped2\n";

				$file1 = $tmp_seq1;
				$file2 = $tmp_seq2;
			}
			else
			{
				my $num_unmapped;

				my $tmp_seq = $list."_read".$j.".gz";

				foreach my $select_read_id (sort keys %select_read)
				{
					# check about read id info
					delete $seq_hash{$select_read_id};	# delete mapped read id : fasta id
				}

				open(TFH, ">:gzip", $tmp_seq) || die "Can not open temp seq (unmapped) file: $tmp_seq\n";
				foreach my $id (sort keys %seq_hash)
				{
					print TFH ">".$id."\n".$seq_hash{$id}."\n";
					$num_unmapped++;
				}
				close(TFH);

				print "No. of unmapped reads: $num_unmapped\n";
				
				$file = $tmp_seq;
			}
		}
		# step 2.5: create header info of SAM file
		else
		{
			my $fha = IO::File->new(">".$head) || die "Can not open header file $head $!\n";
			my $fhb = IO::File->new($output_sam) || die "Can not open output sam file $output_sam $!\n";
			while(<$fhb>)
			{
				if ($_ =~ m/^@/)
				{
					print $fha $_;
				}
				else
				{
					last;
				}
			}
			last;
		}
	}
	#########################################################
	# end cycles for tophat					#
	#########################################################
	
	#########################################################
	# part2	combine result					#
	#########################################################

	# for strand specific library
	if ($split_plus_minus == 1)
	{
		my $output_plus  = $list."_plus.sam";
		my $output_minus = $list."_minus.sam";

		my $pfh = IO::File->new(">".$output_plus)  || die "Can not open plus strand sam file\n";
		my $mfh = IO::File->new(">".$output_minus) || die "Can not open minus strand sam file\n";

		my %plus_strand = (); 
		my %minus_strand = ();

		if ($sequencing_method eq "SS")
		{
			%plus_strand = ('0' => 1, '256' => 1);
			%minus_strand = ('16' => 1, '272' => 1);
		}
		elsif ($sequencing_method eq "PS")
		{
			%plus_strand = (
				'147' => 1, '99' => 1,
                		'145' => 1, '97' => 1,
                		'355' => 1, '403'=> 1,
                		'353' => 1, '401'=> 1);
			%minus_strand = (
				'163' => 1, '83' => 1, 
				'161' => 1, '81' => 1,
				'339' => 1, '419'=> 1,
				'337' => 1, '417'=> 1);
		}
		else
		{
			die "Error at sequencing method $sequencing_method\n";
		}

		for(my $i=0; $i<=$cycle; $i++)
		{
			my $num_mismatch_plus = 0;
			my $num_mismatch_minus = 0;

			my $infile = $list."_m".$i."/mismatch.sam";

			my $ofh3 = IO::File->new($infile)  || die "Can not open output sam file $infile (mismatch $i)\n";
			while(<$ofh3>)
			{
				my @a = split(/\t/, $_);
				unless($_ =~ m/^@/)
				{
					if (defined $minus_strand{$a[1]})
					{
						print $mfh $_; $num_mismatch_minus++;
					}
					elsif (defined $plus_strand{$a[1]}) 
					{
						print $pfh $_; $num_mismatch_plus++;
					} 
					else 
					{
					 	die "Error at line $_";
					}
				}
			}
			$ofh3->close;

			print "Plus strand in mismatch $i-th cycle: $num_mismatch_plus\nMinus strand in mismatch $i-th cycle: $num_mismatch_minus\n";

			my $folder = $list."_m".$i;
			#system("rm -rf $folder");
			print "temp folder $folder is deleted.\n";
		}
		$pfh->close;
		$mfh->close;

		#########################################################
		# sort and add header to file				#
		#########################################################
		my $output_plus_head = $output_plus.".head";
		my $output_minus_head = $output_minus.".head";

		system("cat $head $output_plus > $output_plus_head");
		system("cat $head $output_minus > $output_minus_head");

		my $to1 = $list."_plus_temp.bam";
		my $to2 = $list."_minus_temp.bam";

		my $oo1 = $list."_plus";
		my $oo2 = $list."_minus";

		system("samtools view -bS -o $to1 $output_plus_head");
		system("samtools view -bS -o $to2 $output_minus_head");

		system("samtools sort $to1 $oo1");
		system("samtools sort $to2 $oo2");

		unlink("$output_plus_head");
		unlink("$output_minus_head");
		unlink("$output_plus");
		unlink("$output_minus");
		unlink("$to1");
		unlink("$to2");
		unlink("$head");
	}

	# for unstranded library
	else
	{
		my $output_all = $list."_all.sam";
		my $afh = IO::File->new(">".$output_all) || die "Can not open all sam file\n";

		for(my $i=0; $i<=$cycle; $i++)
		{
			my $aligned_num = 0;
			my $infile = $list."_m".$i."/mismatch.sam";

			my $ofh3 = IO::File->new($infile)  || die "Can not open output sam file $infile (mismatch $i)\n";
			while(<$ofh3>)
			{
				my @a = split(/\t/, $_);
				unless($_ =~ m/^@/)
				{
					print $afh $_; 
					$aligned_num++;
				}
			}
			$ofh3->close;

			print "All mapped read mismatch $i: $aligned_num\n";
                        my $folder = $list."_m".$i;
                        #system("rm -rf $folder");
                        print "temp folder $folder is deleted.\n";
		}
		$afh->close;

		#################################################
		# sort and add header to file			#
		#################################################
		my $output_all_head  = $output_all.".head";
		system("cat $head $output_all > $output_all_head");

		my $to3 = $list."_all_temp.bam";
        	my $oo3 = $list."_all";

		system("samtools view -bS -o $to3 $output_all_head");
		system("samtools sort $to3 $oo3");

		unlink("$head");
		unlink("$output_all");
		unlink("$output_all_head");
		unlink("$to3");
	}
	#########################################################
	# other pipeline 					#
	#########################################################
	print "$list is finished \n";
}

#################################################################
# generate exp_raw info if provide gene position		#
#################################################################
if ($gene_position)
{
	# my $cmd_exp_raw = "get_exp_raw.pl -i $list_file -s $sequencing_method -a $gene_position";
	my $cmd_exp_raw = "${FindBin::RealBin}/get_exp_raw.pl -i $list_file -s $sequencing_method -a $gene_position";
	print "\n\nGet expression raw count ... \n\n$cmd_exp_raw\n";
	system($cmd_exp_raw) && die "Error at cmd $cmd_exp_raw";
	print "\nexpression raw count is generated.\n";

	# my $cmd_uniq_mapped_num = "get_uniq_mapped_num.pl -i $list_file -s $sequencing_method";
	my $cmd_uniq_mapped_num = "${FindBin::RealBin}/get_uniq_mapped_num.pl -i $list_file -s $sequencing_method";
	print "\n\nGet uniq mapped num ... \n\n$cmd_uniq_mapped_num\n";
	system($cmd_uniq_mapped_num) && die "Error at cmd $cmd_uniq_mapped_num";
        print "\nuniq mapped num is generated.\n";

	#########################################################
	# generate exp_rpkm info if provide gene length		#
	#########################################################
	if ($gene_length)
	{
		# generate expression_adjust file;
		my ($libsize, $exp_sense_adjust) = ("", "exp_sense_adjust");
		if (-s "uniq_mapped_num") 
		{
			$libsize = `cut -f2 uniq_mapped_num`;
			chomp($libsize);
			$libsize =~ s/\n/\t/ig;
			$libsize = "libsize\t$libsize";
		}
		else
		{
			print STDERR "Can not get library size from file uniq_mapped_num\n";
			exit;
		}

		my $expraw = IO::File->new("exp_sense_raw") || die "Can not open exp_sense_raw file $!\n";
		my $expadj = IO::File->new(">".$exp_sense_adjust) || die "Can not open exp_sense_adjust file $!\n";
		my $title = <$expraw>;
		print $expadj $title;
		print $expadj $libsize."\n";
		while(<$expraw>) { print $expadj $_; }
		$expadj->close;
		$expraw->close;
			
		# my $cmd_exp_rpkm = "get_exp_rpkm.pl -x $gene_length -e $exp_sense_adjust";
		my $cmd_exp_rpkm = "${FindBin::RealBin}/get_exp_rpkm.pl -x $gene_length -e $exp_sense_adjust";
		print "\n\nGet expression of RPKM ... \n\n$cmd_exp_rpkm\n";
		system($cmd_exp_rpkm) && die "Error at cmd $cmd_exp_rpkm";
		print "\nexpression of RPKM is generated.\n";
	}
}

#################################################################
# kentnf: subroutine						#
#################################################################
=head1 Sub:check_list_file

=cut
sub check_list_file
{
	my ($list_file, $sequencing_method) = @_;
	my %list_read;
	my @file_surfix = (".gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz", ".fq", ".fq.gz", ".fastq", ".fastq.gz");
	my $fh = IO::File->new($list_file) || die "Can not open list file $list_file $!\n";
	my ($file, $file1, $file2);
	while($file = <$fh>)
	{
		chomp($file);

		# check if the read files are exist
		# link the list and read files to hash 
		if ($sequencing_method eq "PS" || $sequencing_method eq "PE")
		{
			$file1 = $file."_1"; 
			$file2 = $file."_2";
			if (-s $file1 && -s $file2)
			{
				$list_read{$file} = $file1."\t".$file2;
			}
			else
			{
				foreach my $surfix (@file_surfix)
				{
					if (-s $file1.$surfix && -s $file2.$surfix)
					{
						$list_read{$file} = $file1.$surfix."\t".$file2.$surfix;
					}
				}
			}
		}
		else
		{
			if (-s $file) 
			{
				$list_read{$file} = $file;
			}
			else
			{
				foreach my $surfix (@file_surfix)
				{
					if (-s $file.$surfix) { 
						$list_read{$file} = $file.$surfix;
					}
				}
			}
		}

		if(defined $list_read{$file}) 
		{
			print STDERR "Locate file for list: $file; File: $list_read{$file}\n";
		}
		else
		{
			print STDERR "Error! Can not locate read for list: $file. (Location: sub->check_list_file) \n";
		}
	}
	return %list_read;
}
=head1 Sub:seq_to_hash

=cut
sub seq_to_hash
{
	my ($seq_file, $cycle) = @_;
	my %seq_hash;
	my ($fh, $format, $id, $seq, $num_of_seq, $read_len, $all_base, $pre_seg_len);
	if ($seq_file =~ m/\.gz/) {
		open($fh, "<:gzip", $seq_file) || die "Location: sub->seq_to_hash. Can not open read file $seq_file $!\n";
	} else {
		open($fh, $seq_file) || die "Location: sub->seq_to_hash. Can not open read file $seq_file $!\n";
	}
	while ($id = <$fh>)
	{
		chomp($id);
		unless($format) 
		{
			if ($id =~ m/^>/) {$format = "fa"; }
			elsif ($id =~ m/^@/) {$format = "fq"; }
			else { die "Location: sub->seq_to_hash. Can not identify the format of read file: $seq_file\n"; }
		}
		if ($format eq "fq") { $id =~ s/^@//; }
		else { $id =~ s/^>//; }
		$id =~ s/\s+.*//;
		$id =~ s/\/\d+//;
		$seq = <$fh>;
		chomp($seq);
		$all_base = $all_base+length($seq);
		if ($id && $seq)
		{
			if ($cycle > 0) { $seq_hash{$id} = $seq;}
			$id = ""; $seq = "";
		}
		else
		{
			print "Location: sub->seq_to_file. Error in read file $seq_file; ID: $id\tSEQ: $seq\n";
		}
		if ($format eq "fq") { <$fh>; <$fh>; }
		$num_of_seq++;
	}
	$fh->close;
	$pre_seg_len = ($all_base/$num_of_seq)/2;
	$pre_seg_len =~ s/\..*//;
	if ($pre_seg_len > 25) { $pre_seg_len = 25; }
	print "Read File: $seq_file; Format: $format; No. of read: $num_of_seq; Predict Segment Length: $pre_seg_len;\n";
	return (\%seq_hash, $pre_seg_len);
}

=head1 Sub:single_end_filter

=cut
sub single_end_filter
{
	my ($insam, $outsam, $mismatch, $cycle) = @_;
	$mismatch =~ s/--read-mismatches //;
	# create the mismatch hash
	my %pattern_mismatch = ();
	for (my $i=0; $i<=$mismatch; $i++)
	{
		$pattern_mismatch{"NM:i:$i"} = 1;
	}

	# proper flag to hash;
	my %flag = ( '16' => 1, '272' => 1, '0' => 1, '256' => 1);

	my %select_read;

	# filter it by flag, and record the mismatch alignment
	my $fout = IO::File->new(">".$outsam) || die "can not open output sam file: $outsam\n";
	my $fhin = IO::File->new($insam) || die "Can not open file $insam\n";
	while(<$fhin>)
	{
		my @a = split(/\t/, $_);

		my $has_pattern = 0; # if has pattern, it is not mismatch.

		# output the temp SAM file base on the flag and mismatch 
		if ($_ =~ m/^@/)
		{
			print $fout $_;
		}
		elsif ( defined $flag{$a[1]}  )
		{
			# find mismatch info base on misatch info
			foreach my $a (@a)
			{
				if (defined $pattern_mismatch{$a})
				{
					$has_pattern = 1;
				}
			}

			if ($has_pattern == 1)
			{
				print $fout $_;
				if ($cycle > 0) {$select_read{$a[0]} = 1;}
			}
		}
		else
		{
			next;
		}
	}
	$fhin->close;
	$fout->close;

	return %select_read;
}

=head1 Sub:paired_end_filter

=cut
sub paired_end_filter
{
	my ($insam, $outsam, $temp, $mismatch, $cycle) = @_;
	$mismatch =~ s/--read-mismatches //;
	# create the mismatch hash
	my %pattern_mismatch = ();
	for (my $i=0; $i<=$mismatch; $i++)
	{
		$pattern_mismatch{"NM:i:$i"} = 1;
	}
	
	# proper flag to hash;
	my %flag = (
		'163' => 1, '83' => 1, '147' => 1, '99' => 1,
		'161' => 1, '81' => 1, '145' => 1, '97' => 1,
		'339' => 1, '419'=> 1, '355' => 1, '403'=> 1,
		'337' => 1, '417'=> 1, '353' => 1, '401'=> 1
	);
	
	# mismatch record hash
	my %mismatch;

	# filter it by flag, and record the mismatch alignment
	my $ftmp = IO::File->new(">".$temp) || die "can not open temp file: $temp\n";
	my $fhin = IO::File->new($insam) || die "Can not open file $insam\n";
	while(<$fhin>)
	{
		my @a = split(/\t/, $_);

		my $has_pattern = 0; # if has pattern, it is not mismatch.

		# output the temp SAM file base on the flag 
		if ($_ =~ m/^@/)
		{
			print $ftmp $_;
		}
		elsif ( defined $flag{$a[1]} )
		{
			print $ftmp $_;

			# find mismatch info base on misatch info
			foreach my $a (@a)
			{
				if (defined $pattern_mismatch{$a})
				{
					$has_pattern = 1;
				}
			}

			if ($has_pattern == 0)
			{
				$mismatch{$a[0]."#".$a[2]."#".$a[3]."#".$a[7]} = 1;
				$mismatch{$a[0]."#".$a[2]."#".$a[7]."#".$a[3]} = 1;
			}
		}
		else
		{
			next;
		}
	}
	$fhin->close;
	$ftmp->close;

	my %select_read;
	# filter the mismatch alignment base on mismatch record
	my $ftmp = IO::File->new($temp) || die "can not open temp file: $temp\n";	
	my $fout = IO::File->new(">".$outsam) || die "Can not open output file: $outsam\n";
	while(<$ftmp>)
	{
		chomp;
		my @a = split(/\t/, $_);

		if ($_ =~ m/^@/)
		{
			print $fout $_."\n";
		}		
		elsif ( defined $mismatch{$a[0]."#".$a[2]."#".$a[3]."#".$a[7]} || defined $mismatch{$a[0]."#".$a[2]."#".$a[7]."#".$a[3]} )
		{
			next;
		}
		else
		{
			if ($cycle > 0) { $select_read{$a[0]} = 1; }
			print $fout $_."\n";
		}
	}
	$ftmp->close;
	$fout->close;

	return %select_read;
}

=head1 Sub:parse_cigar

=cut
sub parse_cigar
{
	my $cigar = shift;

	my $str_len = length($cigar);

	my $num = "";; my $mapped_length = 0;

	for(my $i=0; $i<$str_len; $i++)
	{
		my $str = substr($cigar, $i, 1);

		if ($str =~ m/\d+/)
		{
			$num = $num.$str;
		}
		elsif ($str eq "M" || $str eq "N" || $str eq "I")
		{
			$mapped_length = $mapped_length + $num;
			$num = "";
			
		}
		elsif ($str eq "D") 
		{
			$num = "";
		}
	}

	return $mapped_length;
}
