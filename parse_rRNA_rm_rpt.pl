#!/usr/bin/perl

use IO::File;

my $usage = qq'
perl parse_rRNA_rm_rpt.pl  rRNA_rm.report

';

my $file = shift || die $usage;

print "#sample\ttotal\trRNA\tclean\n";

my $fh = IO::File->new($file) || die "Can not open rRNA_rm report file $file $!\n";
while(<$fh>)
{
	chomp;
	#if ($_ =~ m/^Locate/ || $_ =~ m/^Read/) {  $start = 0; next; }
	if ( $_ =~ m/^bowtie/ )  { $start = 1; }
	else { $start = 0; next; }

	if ($start == 1)
	{
		$sample = $_;			@a = split(/\s+/, $sample);
		$total = <$fh>;  chomp($total);	@b = split(/\s+/, $total);
		$align = <$fh>;  chomp($align);	@c = split(/\s+/, $align);
		$failed = <$fh>; chomp($failed);@d = split(/\s+/, $failed);
		<$fh>;
		print $a[8]."\t".$b[3]."\t".$c[8]."\t".$d[6]."\n";
	}
}
$fh->close;
