#!/usr/bin/perl

use strict;
use warnings;
use IO::File;

# input is the raw count, rpkm file, pairwise comparison list, and the sample name must be sample_replicateN;

my $usage = qq'
Perl TSlimma_pipeline.pl  raw_count  rpkm_file  comparison_list  output

';

my $raw_count = shift || die $usage;
my $rpkm_file = shift || die $usage;
my $comp_list = shift || die $usage;
my $output = shift || 'DE_out';
#================================================================
# store comparison list to hash	
# key: compN; value: sample1 vs sample2				
# file format:
# Time1_sample \t Time2_sample \t ... \t TimeN_sample \n; 
#================================================================
my @comparison;

my $fh = IO::File->new($comp_list) || die "Can not open comparison list file $comp_list $!\n";
while(<$fh>)
{
	chomp;
	push(@comparison, $_);
}
$fh->close;

#================================================================
# parse raw count dataset					
# save comparison data files without zero
#================================================================
my %raw; my %title;

my $fh = IO::File->new($raw_count) || die "Can not open raw count file $raw_count $!\n";
my $title = <$fh>;
chomp($title);
my @t = split(/\t/, $title);
for(my $i=1; $i<@t; $i++)
{
	my $sample_name = $t[$i];
	my $sample = $sample_name;
	$sample =~ s/_rep\d+//;
	if (defined $title{$sample} ) {
		$title{$sample}.="\t".$sample_name;
	} else {
		$title{$sample} = $sample_name;
	}
}

while(<$fh>)
{
	chomp;
	my @a = split(/\t/, $_);
	my $gene = $a[0];
	for(my $i=1; $i<@a; $i++)
	{
		$sample_name = $t[$i];
		$sample = $sample_name;
		$sample =~ s/_rep\d+//;
		my $raw_value = $a[$i];
		if (defined $raw{$gene}{$sample}) {
			$raw{$gene}{$sample} = $raw{$gene}{$sample}."\t".$raw_value;
		} else {
			$raw{$gene}{$sample} = $raw_value;
		}
	}
}
$fh->close;

#================================================================
# parse RPKM file
# get mean and ratio from RPKM
#================================================================
my %total;      # the total expression of all rep for one sample
my %RPKM;       # the rpkm
my $fh = IO::File->new($rpkm_file) || die "Can not open RPKM file $rpkm_file $!\n";
my $title_R = <$fh>;
chomp($title_R);
my @tr = split(/\t/, $title_R);
while(<$fh>)
{
        chomp;
        my @a = split(/\t/, $_);
        my $gene = $a[0];
        for(my $i=1; $i<@a; $i++)
        {
                $sample_name = $tr[$i]; # sample name
                $sample = $sample_name;
                $sample =~ s/_rep\d+//;
                $RPKM_value = $a[$i];
                if (defined $total{$gene}{$sample})
                {
                        $total{$gene}{$sample} = $total{$gene}{$sample} + $RPKM_value;
                }
                else
                {
                        $total{$gene}{$sample} = $RPKM_value;
                }


                if (defined $RPKM{$gene}{$sample})
                {
                        $RPKM{$gene}{$sample}.="\t".$RPKM_value;
                }
                else
                {
                        $RPKM{$gene}{$sample} = $RPKM_value;
                }
        }
}
$fh->close;

# get adjusted P value 
my %padj; my %replicate; my %comp_sample;

foreach my $comp (@comparison)
{
	my @sample = split(/\t/, $comp);

	$time++;
	my $raw_file = "TimeSeries".$time."_raw"; 
	my $vsd_file = "TimeSeries".$time."_vsd";
	my $out_file = "TimeSeries".$time."_out";
	my $deg_file = "TimeSeries".$time."_DEG";
	my $fit_file = "TimeSeries".$time."_DEG_Filter";

	my @v;

	# generate raw count file for TimeSeries comparison and remove gene without any count
	my $rfh = IO::File->new(">".$raw_file) || die "Can not open file $raw_file\n";

	# print title
	print $rfh "gene";
	foreach my $s (@sample) 
	{ 
		print $rfh "\t".$title{$s};
		my @t = split(/\t/, $title{$s});
		$replicate{$s} = scalar(@t);
	}
	print $rfh "\n";

	# print raw count
	foreach my $gene (sort keys %raw) 
	{
		my $sum = 0;
		my $line = $gene;

		foreach my $s (@sample)
		{
			@v = split(/\t/, $raw{$gene}{$s});
			foreach my $v1 (@v) { $sum = $sum + $v1; }
			$line.="\t".$raw{$gene}{$s};
		}
		$line.="\n";

		if ($sum > 0) { print $rfh $line; }
	}
	$rfh->close;

	# VST transform
	my $r1 = vsd_transform($raw_file, $vsd_file, \@sample, \%replicate);
	my $tmp1 = IO::File->new(">temp1.R") || die "Can not open temp1.R file $!\n";
        print $tmp1 $r1;
        $tmp1->close;
        system("R --no-save < temp1.R") && die "Error at cmd R --no-save < temp1.R\n";

	# add gene title and remove " from vsd file
	fix_vsd($vsd_file);

	# DEG analysis using limma
	my $r2 = limma_code($vsd_file, $out_file, \@sample, \%replicate);
	my $tmp2 = IO::File->new(">temp2.R") || die "Can not open temp2.R file $!\n";
	print $tmp2 $r2;
	$tmp2->close;
	system("R --no-save < temp2.R") && die "Error at cmd R --no-save < temp2.R\n";

	# parse R output file and save 1) F value, 2) P value, and 3) adjusted p value to hash
	my $ofh = IO::File->new($out_file) || die "Can not open DESeq output file $out_file $!\n";
	<$ofh>;
	while(<$ofh>)
	{
		chomp;
		$_ =~ s/"//ig;
		my @a = split(/\t/, $_);
		my $last_array = scalar(@a);
		($gid, $f_value, $p_value, $adj_p_value) = ($a[1], $a[$last_array-3], $a[$last_array-2], $a[$last_array-1]);
		$padj{$gid}{$comp} = $f_value."\t".$p_value."\t".$adj_p_value;
	}
	$ofh->close;

	# output RPKM, mean of RPKM, and padj hash (include F value, P value, and adjusted P value)
	my $dfh = IO::File->new(">".$deg_file) || die "Can not open DEG file $deg_file $!\n";
	my $ffh = IO::File->new(">".$fit_file) || die "Can not open DEG filter file $fit_file $!\n"; 

	print $dfh "GeneID";
	print $ffh "GeneID";
	foreach my $s (@sample)
	{
		print $dfh "\t".$title{$s}."\tmean";
		print $ffh "\t".$title{$s}."\tmean";
	}
	print $dfh "\tF\tP.Value\tadj.P.Val\n";
	print $ffh "\tF\tP.Value\tadj.P.Val\n";

	my $out_line; my @b; my $num_of_deg = 0;
	foreach my $gene (sort keys %RPKM)
	{
		@b = ();
		$out_line = $gene;
		foreach my $s (@sample)
		{
			my $mean = $total{$gene}{$s}/$replicate{$s};
			$mean = sprintf("%.2f", $mean);
			$out_line.="\t".$RPKM{$gene}{$s}."\t".$mean;
		}
		
		if (defined $padj{$gene}{$comp})
		{
			$out_line.="\t".$padj{$gene}{$comp}."\n";
			@b = split(/\t/, $padj{$gene}{$comp});
		}
		else
		{
			$out_line.="\tNA\tNA\tNA\n";
			
		}

		if ($b[2] > 0 && $b[2] < 0.05)
		{
			print $ffh $out_line;
			$num_of_deg++;
		}
		print $dfh $out_line;
	}

	$dfh->close;
	$ffh->close;

	print $comp."\t".$num_of_deg."\n";
}

die;

my %ratio;
foreach my $gene (sort keys %mean)
{
	foreach my $comparison (@comparison)
	{
		my ($compA, $compB) =  split(/\t/, $comparison);
		$meanA = $mean{$gene}{$compA};
		$meanB = $mean{$gene}{$compB};
		if ($low_RPKM)
		{
			if ($compA == 0 && $compB == 0 )
			{
				$ratio = 1;
			}
			elsif ($compA == 0)
			{
				$ratio = $compB / $low_RPKM;
			}
			elsif ($compB == 0)
			{
				$ratio = $low_RPKM / $compA;
			}
			else
			{
				$ratio = $compB / $compA;
			}
		}
		else
		{
			if ($meanA == 0 && $meanB == 0 )
			{
				$ratio = 1;
			}
			elsif ($meanA == 0)
			{
				$ratio = $meanB / 0.01;
			}
			elsif ($meanB == 0)
			{
				$ratio = 0.01 / $meanA;
			}
			else
			{
				$ratio = $meanB / $meanA;
			}
		}
		$ratio{$gene}{$comparison} = $ratio;
	}
}

#================================================================
# out put result
#================================================================
my $out1 = IO::File->new(">".$output) || die "Can not open output file $output\n";
my $out2 = IO::File->new(">".$output."_filter") || die "Can not open filtered output file \n";


print $out1 "GeneID";
print $out2 "GeneID";
foreach my $comp (@comparison)
{
	my ($sampleA, $sampleB) = split(/\t/, $comp);
	print $out1 "\t".$title{$sampleA}."\tmean\t".$title{$sampleB}."\tmean\tratio\tadjust p";
	print $out2 "\t".$title{$sampleA}."\tmean\t".$title{$sampleB}."\tmean\tratio\tadjust p";
}
print $out1 "\n";
print $out2 "\n";

my ($out_line, $sig); 
my %report;
# set these parameters before
$fc1=2; $fc2 = 0.5; $padj_cutoff = 0.05;

foreach my $gene (sort keys %RPKM)
{
	$out_line = $gene;
	$sig = 0;
	
	foreach my $comp ( @comparison )
	{
		my ($sampleA, $sampleB) = split(/\t/, $comp);
	
	      	#$out_line.="\t".$RPKM{$gene}{$sampleA}."\t".
		#	$mean{$gene}{$sampleA}."\t".
		#	$RPKM{$gene}{$sampleB}."\t".
		#	$mean{$gene}{$sampleB}."\t".
		#	$ratio{$gene}{$comp}."\t";

		my $meanA = $mean{$gene}{$sampleA};
		$meanA = sprintf("%.2f", $meanA);
		my $meanB = $mean{$gene}{$sampleB};
                $meanB = sprintf("%.2f", $meanB);
		my $ratio = $ratio{$gene}{$comp};
		$ratio = sprintf("%.2f", $ratio);

		$out_line.="\t".$RPKM{$gene}{$sampleA}."\t".
                        $meanA."\t".
                        $RPKM{$gene}{$sampleB}."\t".
                        $meanB."\t".
                        $ratio."\t";

		if (defined $padj{$gene}{$comp}) {
			$out_line.=$padj{$gene}{$comp};
		} else {
			$out_line.="NA";
		}

		if (defined $padj{$gene}{$comp}) {
			if (($ratio{$gene}{$comp} > $fc1 || $ratio{$gene}{$comp} < $fc2) && $padj{$gene}{$comp} < $padj_cutoff)
			{
				$sig = 1;
				$report{$comp}++;
			}
		}
	}
	if ( $sig ) { print $out2 $out_line."\n"; }
	print $out1 $out_line."\n";
}
$out1->close;
$out2->close;

# report the number of DE genes for every comparison
foreach my $comp (sort keys %report)
{
	$num = $report{$comp};
	$comp =~ s/\s+/ vs /;
	print $comp."\t".$num."\n";
}

#================================================================
# kentnf: subroutine
#================================================================
sub vsd_transform
{
	my ($input, $output, $sample, $num ) = @_;
	my $pwd = `pwd`;
	chomp($pwd);

	my $factor;
	foreach my $sample_name ( @$sample )
	{
		my $num_of_replicate = $$num{$sample_name};
		for(my $i=0; $i<$num_of_replicate; $i++) { $factor.=" \"$sample_name\","; }
	}
	$factor =~ s/,$//;

	my $r_code = qq'
setwd(\'$pwd\')
library(DESeq)
countsTable<-read.delim("$input", header=TRUE, stringsAsFactors=TRUE)
rownames(countsTable)<-countsTable\$gene
countsTable<-countsTable[, -1]
conds <- factor( c( $factor ) )
cds<-newCountDataSet(countsTable, conds)
cds <- estimateSizeFactors( cds )
sizeFactors( cds )
cdsBlind <- estimateDispersions( cds, method="blind" )
vsd <- getVarianceStabilizedData( cdsBlind )
write.table(vsd, sep="\\t", file="$output")
';
	return $r_code;	
}

sub limma_code
{
	my ($input, $output, $sample, $num ) = @_;
	my $pwd = `pwd`;
	chomp($pwd);

	my $factor;
	my $design;
	my $comparison;

	my @sample = @$sample;
	for(my $i=0; $i<@sample; $i++)
	{
		my $sample_name = $sample[$i];
		
		if ($i > 0)
		{
			my $comp = $sample[$i]."-".$sample[$i-1];
			$comparison.= "\"$comp\"".", ";
		}

		$design.="\"$sample_name\", ";

                my $num_of_replicate = $$num{$sample_name};

                for(my $j=0; $j<$num_of_replicate; $j++) { my $k = $i+1; $factor.= $k.","; }
        }

	$design =~ s/, $//;
        $factor =~ s/,$//;

	my $xcomp;

	my %uniq;

	foreach my $aa (@sample)
	{
		foreach my $bb (@sample)
		{
			if ($aa ne $bb)
			{
				unless($uniq{$bb.$aa})
				{
					$xcomp.= "\"".$bb."-".$aa."\", ";
					$uniq{$aa.$bb} = 1;
					$uniq{$bb.$aa} = 1;
				}
			}
		}
	}

	my $r_code = qq'
setwd(\'$pwd\')
library(limma)

eset <- read.table(file="$input", header=TRUE)
rownames(eset)<-eset\$gene
eset<-eset[,-1]

design <- model.matrix(~ -1+factor(c($factor)))
colnames(design) <- c($design)
fit <- lmFit(eset, design)

# Which genes respond (i.e., change over time) ?
cont.f <- makeContrasts($xcomp levels=design)
fitTS <- contrasts.fit(fit, cont.f)
fitTS <- eBayes(fitTS)
compTimeF <- topTableF(fitTS, adjust="BH", number=50000)
write.table(compTimeF, sep="\\t", file="$output")
';
	return $r_code;

}

sub fix_vsd
{
	my $infile = shift;
	my $out = IO::File->new(">temp.vsd") || die "Can not open temp vsd file $!\n";
	my $in = IO::File->new($infile) || die "Can not open vsd file $infile $!\n";
	my $title = <$in>;
	$title =~ s/"//ig;
	print $out "gene\t".$title;
	while(<$in>)
	{
		$_ =~ s/"//ig;
		print $out $_;
	}
	$in->close;
	$out->close;
	system("mv temp.vsd $infile") && die "Error in command : mv temp.vsd $infile $!\n";
}
