#!/usr/bin/perl

my $usage = qq'
USAGE: perl $0  list  position_file

';

my $list = shift || die $usage;
my $position = shift || die $usage;

my $s = "SS";
my $position_line = 40000;

my $line = `wc -l $position`; chomp $line;
my $cycle = int( $line / $position_line ) + 1;

print "position line:$line\tcycle: $cycle\n";

for(my $i=1; $i<=$cycle; $i++)
{
	my $pos = $position.$i;
	
	generate_pos($position, $i, $position_line, $pos);

	my $out = "exp$i";
	my $cmd = "get_exp_raw.pl -i $list -s $s -a $pos -o $out";
	print $cmd."\n";
	system($cmd) && die "Error in command $cmd \n";
	unlink($pos);
}

sub generate_pos
{
	my ($position, $i, $position_line, $pos) = @_;

	my $n = 0;
	my $end = $i * $position_line;
	my $start = $end - $position_line + 1;	
	open(OUT, ">".$pos) || die $!;
	open(FH, $position) || die $!;
	while(<FH>)
	{
		chomp;
		$n++;
		if ($n>=$start && $n <= $end) { print OUT $_."\n"; }
	}
	close(FH);
	close(OUT);
}



#Perl get_exp_raw.pl -i read_list -s SS -a tomato_gene_position



