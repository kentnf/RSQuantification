#!/usr/bin/perl

=head

05/23/2013
combine DESeq and edgeR to one program

06/12/2012
generate pdf for each comparison, including PCA plot, MA plot, and BCV plot.

06/11/2012
identify DE genes using edgeR

=cut

use strict;
use warnings;
use IO::File;
use Getopt::Long;

my $usage = qq'
Perl DEGs_pipeline.pl [options]

	-i raw_count 
	-r rpkm_file 
	-c comparison_list 
	-p program (default = DESeq)
	-o output (default = program)

* the sample name for raw count and rpkm file should be like sample_repN
* the program should be DESeq or edgeR

';

my ($help, $raw_count, $rpkm_file, $comp_list, $program, $output);

GetOptions(
	"h"	=> \$help,
	"i=s"	=> \$raw_count,
	"r=s"	=> \$rpkm_file,
	"c=s"	=> \$comp_list,
	"p=s"	=> \$program,
	"o=s"	=> \$output
);

die $usage if $help;
die $usage unless $raw_count;
die $usage unless $rpkm_file;
die $usage unless $comp_list;
$program ||= 'DESeq';
$output ||= $program;

if ($program eq 'DESeq' || $program eq 'edgeR') {} else { die "Error in program\n"; }

#================================================================
# store comparison list to hash	
# key: compN; value: sample1 vs sample2				
# file format:
# sample1 \t sample2 \n; 
#================================================================
my @comparison = comparison_to_array($comp_list);

sub comparison_to_array
{
	my $comp_list = shift;
	my @comparison;
	my $fh = IO::File->new($comp_list) || die "Can not open comparison list file $comp_list $!\n";
	while(<$fh>) { chomp; push(@comparison, $_); }
	$fh->close;
	return @comparison;
}

#================================================================
# parse raw count dataset					
# save comparison data files without zero
#================================================================
my ($title, $raw) = raw_count_to_hash($raw_count);

sub raw_count_to_hash
{
	my $raw_count = shift;

	my %raw; my %title; my ($gene, $sample_name, $sample, $raw_value);

	my $fh = IO::File->new($raw_count) || die "Can not open raw count file $raw_count $!\n";

	# parse raw count title
	my $title = <$fh>; chomp($title);
	my @t = split(/\t/, $title);
	for(my $i=1; $i<@t; $i++)
	{
		$sample_name = $t[$i];
		$sample = $sample_name;
		$sample =~ s/_rep\d+//;
		if (defined $title{$sample} ) {
			$title{$sample}.="\t".$sample_name;
		} else {
			$title{$sample} = $sample_name;
		}
	}

	# parse raw count value
	while(<$fh>)
	{
		chomp;
		my @a = split(/\t/, $_);
		$gene = $a[0];
		for(my $i=1; $i<@a; $i++)
		{
			$sample_name = $t[$i];
			$sample = $sample_name;
			$sample =~ s/_rep\d+//;
			$raw_value = $a[$i];
			if (defined $raw{$gene}{$sample}) {
				$raw{$gene}{$sample} = $raw{$gene}{$sample}."\t".$raw_value;
			} else {
				$raw{$gene}{$sample} = $raw_value;
			}
		}
	}
	$fh->close;

	return(\%title, \%raw);
}

#================================================================
# save comparison data files without zero
# statistics analysis 
#================================================================
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
	my @v1; my @v2;
	my $rfh = IO::File->new(">".$raw_file) || die "Can not open file $raw_file\n";
	print $rfh "gene\t",$$title{$sampleA},"\t",$$title{$sampleB},"\n";
	foreach my $gene (sort keys %$raw) {
		# remove zero
		my $sum = 0;
		@v1 = split(/\t/, $$raw{$gene}{$sampleA});
		@v2 = split(/\t/, $$raw{$gene}{$sampleB});
		$replicate{$sampleA} = scalar(@v1);
		$replicate{$sampleB} = scalar(@v2);
		foreach my $v1 (@v1) { $sum = $sum + $v1; }
		foreach my $v2 (@v2) { $sum = $sum + $v2; }
		if ($sum > 0) {
			print $rfh $gene,"\t",$$raw{$gene}{$sampleA},"\t",$$raw{$gene}{$sampleB}."\n";
		}
	}
	$rfh->close;

	my $r; my $gene_column; my $pvalue_column; 

	if ($program eq 'DESeq')
	{
		$r = generate_r_deseq($raw_file, $out_file, $sampleA, $sampleB, scalar(@v1), scalar(@v2));
		$pvalue_column = 8;
		$gene_column = 1;
	}
	elsif ($program eq 'edgeR')
	{
		$r = generate_r_edger($raw_file, $out_file, $sampleA, $sampleB, scalar(@v1), scalar(@v2));
		$pvalue_column = 4;
		$gene_column = 0;
	}

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
		$padj{$a[$gene_column]}{$comp} = $a[$pvalue_column];
	}
	$ofh->close;
}

#================================================================
# parse RPKM file
# get mean and ratio from RPKM
#================================================================
# total : key: gene_id sample, value: the total expression of all rep for one sample
# rpkm  : key: gene_id sample, value: RPKM for each replictes
my ($total, $RPKM) = rpkm_to_hash($rpkm_file);

sub rpkm_to_hash
{
	my $rpkm_file = shift;

	my $fh = IO::File->new($rpkm_file) || die "Can not open RPKM file $rpkm_file $!\n";

	my %RPKM; my %total; my ($gene, $sample_name, $sample, $RPKM_value);

	# parse title
	my $title_R = <$fh>; chomp($title_R);
	my @tr = split(/\t/, $title_R);

	# parse rpkm value
	while(<$fh>)
	{
		chomp;
		my @a = split(/\t/, $_);
		$gene = $a[0];
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

	return(\%total, \%RPKM);
}

my %mean;
foreach my $gene (sort keys %$total)
{
	foreach my $sample ( sort keys $$total{$gene} )
	{
		if (defined $comp_sample{$sample})
		{
			my $total = $$total{$gene}{$sample};
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

my %ratio;  my $low_RPKM;
foreach my $gene (sort keys %mean)
{
	foreach my $comparison (@comparison)
	{
		my ($compA, $compB) =  split(/\t/, $comparison);
		my $meanA = $mean{$gene}{$compA};
		my $meanB = $mean{$gene}{$compB};
		my $ratio;
		if (defined $low_RPKM && $low_RPKM > 0 )
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
	print $out1 "\t".$$title{$sampleA}."\tmean\t".$$title{$sampleB}."\tmean\tratio\tadjust p";
	print $out2 "\t".$$title{$sampleA}."\tmean\t".$$title{$sampleB}."\tmean\tratio\tadjust p";
}
print $out1 "\n";
print $out2 "\n";

my ($out_line, $sig); 
my %report;
# set these parameters before
my $fc1=2; my $fc2 = 0.5; my $padj_cutoff = 0.05;

foreach my $gene (sort keys %$RPKM)
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

		$out_line.="\t".$$RPKM{$gene}{$sampleA}."\t".
                        $meanA."\t".
                        $$RPKM{$gene}{$sampleB}."\t".
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
	my $num = $report{$comp};
	$comp =~ s/\s+/ vs /;
	print $comp."\t".$num."\n";
}

#================================================================
# kentnf: subroutine
#================================================================
sub generate_r_deseq
{
	my ($input, $output, $sampleA, $sampleB, $numA, $numB ) = @_;
	my $pwd = `pwd`;
	chomp($pwd);

	my $DispEsts_pdf = $sampleA."_".$sampleB."_DispEsts.pdf";
	my $DE_pdf = $sampleA."_".$sampleB."_DE.pdf";
	my $hist_pdf = $sampleA."_".$sampleB."_hist.pdf";

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

#plot DispEsts pdf
plotDispEsts <- function( cds )
{
plot(
rowMeans( counts( cds, normalized=TRUE ) ),
fitInfo(cds)\$perGeneDispEsts,
pch = \'.\', log="xy" )
xg <- 10^seq( -.5, 5, length.out=300 )
lines( xg, fitInfo(cds)\$dispFun( xg ), col="red" )
}
pdf("$DispEsts_pdf", width=8, height=6)
plotDispEsts( cds )

comp <- nbinomTest( cds, "$sampleA", "$sampleB" )
write.table( comp, sep="\\t", file="$output" )

#plot DE pdf
plotDE <- function( comp )
plot(
comp\$baseMean,
comp\$log2FoldChange,
log="x", pch=20, cex=.3,
col = ifelse( comp\$padj < .05, "red", "black" ) )

pdf("$DE_pdf", width=8, height=6)
plotDE( comp )

#plot hist pdf
pdf("$hist_pdf", width=8, height=6)
hist(comp\$pval, breaks=100, col="skyblue", border="slateblue", main="")
';
	return $r_code;	
}

sub generate_r_edger
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
detags <- rownames(topTags(et, n = 40000)\$table)
pdf("$ma_pdf",width=8,height=6)
plotSmear(et, de.tags=detags)
abline(h = c(-2, 2), col = "dodgerblue")
';
	return $r_code;
}

