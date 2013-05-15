#!/usr/bin/perl

=head

06/12/2012
generate pdf for each comparison, including PCA plot, MA plot, and BCV plot.

06/11/2012 
identify DE genes using edgeR

=cut

use strict;
use warnings;
use IO::File;

# input is the raw count, rpkm file, pairwise comparison list, and the sample name must be sample_replicateN;
my $usage = qq'
Perl edgeR_pipeline.pl  raw_count  rpkm_file  comparison_list  output

';

my $raw_count = shift || die $usage;
my $rpkm_file = shift || die $usage;
my $comp_list = shift || die $usage;
my $output = shift || 'DE_out';
#================================================================
# store comparison list to hash	
# key: compN; value: sample1 vs sample2				
# file format:
# sample1 \t sample2 \n; 
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

my %padj; my %replicate; my %comp_sample;

foreach my $comp (@comparison)
{
	my ($sampleA, $sampleB) = split(/\t/, $comp);
	$comp_sample{$sampleA} = 1;
	$comp_sample{$sampleB} = 1;
	my $raw_file = $comp; my $out_file = $comp;
	$raw_file =~ s/\t/_/; $out_file =~ s/\t/_/;
	$raw_file = $raw_file."_raw";
	$out_file = $out_file."_out";
	my @v1; @v2;
	my $rfh = IO::File->new(">".$raw_file) || die "Can not open file $raw_file\n";
	print $rfh "gene\t",$title{$sampleA},"\t",$title{$sampleB},"\n";
	foreach my $gene (sort keys %raw) {
		# remove zero
		my $sum = 0;
		@v1 = split(/\t/, $raw{$gene}{$sampleA});
		@v2 = split(/\t/, $raw{$gene}{$sampleB});
		$replicate{$sampleA} = scalar(@v1);
		$replicate{$sampleB} = scalar(@v2);
		foreach my $v1 (@v1) { $sum = $sum + $v1; }
		foreach my $v2 (@v2) { $sum = $sum + $v2; }
		if ($sum > 0) {
			print $rfh $gene,"\t",$raw{$gene}{$sampleA},"\t",$raw{$gene}{$sampleB}."\n";
		}
	}
	$rfh->close;

	my $r = generate_r($raw_file, $out_file, $sampleA, $sampleB, scalar(@v1), scalar(@v2));
	my $tmp = IO::File->new(">temp.R") || die "Can not open temp.R file $!\n";
	print $tmp $r;
	$tmp->close;
	system("R --no-save < temp.R") && die "Error at cmd R --no-save < temp.R\n";

	# parse R output file and save adjusted p value to hash
	my $ofh = IO::File->new($out_file) || die "Can not open DESeq output file $out_file $!\n";
	<$ofh>;
	while(<$ofh>)
	{
		chomp;
		$_ =~ s/"//ig;
		my @a = split(/\t/, $_);
		$padj{$a[0]}{$comp} = $a[4];
	}
	$ofh->close;
}

#================================================================
# parse RPKM file
# get mean and ratio from RPKM
#================================================================
my %total;	# the total expression of all rep for one sample
my %RPKM;	# the rpkm
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

my %mean;
foreach my $gene (sort keys %total)
{
	foreach my $sample (sort keys %{$total{$gene}})
	{
		if (defined $comp_sample{$sample})
		{
			my $total = $total{$gene}{$sample};
			my $rep = $replicate{$sample};

			if ($rep == 0) {
				die "Error at replicate num of $sample\n";
			}

			my $mean;
			if ($total > 0) { $mean = $total/$rep; }
			else		{ $mean = 0; }
			$mean{$gene}{$sample} = $mean;
		}
	}
}

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
			if (($ratio{$gene}{$comp} > $fc1 || $ratio{$gene}{$comp} < $fc2) && $padj{$gene}{$comp} < $padj_cutoff) {
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
=head
sub generate_r
{
	my ($input, $output, $sampleA, $sampleB, $numA, $numB ) = @_;
	my $pwd = `pwd`;
	chomp($pwd);	
	my ($factorA, $factorB);
	for(my $i=0; $i<$numA; $i++) { $factorA.=" \"$sampleA\","; }
	for(my $i=0; $i<$numB; $i++) { $factorB.=" \"$sampleB\","; }
	$factorB =~ s/,$//;
	my $r_code = qq'
setwd(\'$pwd\')
library(DESeq)
countsTable<-read.delim("$input", header=TRUE, stringsAsFactors=TRUE)
rownames(countsTable)<-countsTable\$gene
countsTable<-countsTable[, -1]
conds <- factor( c( $factorA $factorB ) )
cds<-newCountDataSet(countsTable, conds)
cds <- estimateSizeFactors( cds )
cds <- estimateDispersions( cds )
comp <- nbinomTest( cds, "$sampleA", "$sampleB" )
write.table( comp, sep="\\t", file="$output" )
';
	return $r_code;	
}
=cut

sub generate_r
{
	my ($input, $output, $sampleA, $sampleB, $numA, $numB ) = @_;
	my $num_end = 2+$numA+$numB-1; 
        my $pwd = `pwd`;
        chomp($pwd);    
 	
	my $pca_pdf = $sampleA."_".$sampleB."_PCA.pdf";
	my $bcv_pdf = $sampleA."_".$sampleB."_BCV.pdf";
	my $ma_pdf = $sampleA."_".$sampleB."_MA.pdf";

	my ($factorA, $factorB);
        for(my $i=0; $i<$numA; $i++) { $factorA.=" \"$sampleA\","; }
        for(my $i=0; $i<$numB; $i++) { $factorB.=" \"$sampleB\","; }
        $factorB =~ s/,$//;
        my $r_code = qq'
setwd(\'$pwd\')
library(edgeR)
library(limma)
raw.data <- read.delim("$input")
#names(raw.data)

# normalization and filtering
d <- raw.data[, 2:$num_end]
rownames(d) <- raw.data[, 1]
group <- c(rep("$sampleA", $numA), rep("$sampleB", $numB))
d <- DGEList(counts = d, group = group)
dim(d)
cpm.d <- cpm(d)
d <- d[ rowSums(cpm.d > 1) >=3, ]
d <- calcNormFactors(d)

# Data exploration, generate PCA pdf
pdf("$pca_pdf",width=8,height=6)
plotMDS(d, xlim=c(-1,1), labels = c( $factorA $factorB ))

# Estimating the dispersion
d <- estimateCommonDisp(d, verbose=TRUE)
d <- estimateTagwiseDisp(d)
pdf("$bcv_pdf",width=8,height=6)
plotBCV(d)
et <- exactTest(d)
result <- topTags(et, n=40000, adjust.method="BH", sort.by="p.value")
write.table( result, sep="\\t", file="$output" )

# generate MA (Smear) plot
detags <- rownames(topTags(et, n = 40000)$table)
pdf("$ma_pdf",width=8,height=6)
plotSmear(et, de.tags=detags)
abline(h = c(-2, 2), col = "dodgerblue")
';
        return $r_code;
}
