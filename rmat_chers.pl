use warnings;
use strict;
use Getopt::Long;

# Karsten Hokamp, Trinity College Dublin, 2015
# Perl script that detects CHip-Enriched Regions (CHERs) from RMAT output
# input consists of BED files containing chr, start and end position of a probe,
# as well as the associated MAT score (see example lines below)
# output will be written to files with the same name as input files but .cher instead of .bed ending
# values from overlapping probes are averaged using the mean
# a sliding window approach is used for averging scores
# contiguous regions with score above threshold will be reported as CHERs

my $window_size = 300;
my $min_width = 7 * 4 + 24;
my $min_score = 1.5;
my $user_chr = '';

&GetOptions(
    'window_size=i' => \$window_size,
    'min_score=f' => \$min_score,
    'min_width=i' => \$min_width,
    'chr=s' => \$user_chr,
    );

foreach my $file (@ARGV) {

    print STDERR "Processing $file on ".(localtime)."\n";

    open (IN, $file)
	or die;

#chrI    0       24      -6.72839472572012
#chrI    4       28      -6.80581739601674
#chrI    8       32      -6.91517933147194

    # read in all scores 
    my %in = ();
    while (<IN>) {
	chomp;
	my ($chr, $start, $end, $score) = split /\t/, $_;
	if ($user_chr ne '') {
	    next unless ($chr eq $user_chr);
	}
			      
	$start++;
	$end++;
	
	foreach ($start..$end) {
	    push @{$in{$chr}{$_}}, $score;
	}
    }
    close IN;

    my $out_file = $file;
    $out_file =~ s/\.bed$//;
    $out_file .= '.cher';
    open (OUT, ">$out_file")
	or die;
    select OUT;

    my @chr = sort keys %in;
    
    # average overlapping probes
    foreach my $chr (@chr) {
	foreach my $pos (sort { $a <=> $b } keys %{$in{$chr}}) {
	    my $mean = &mean(@{$in{$chr}{$pos}});
	    $in{$chr}{$pos} = $mean;
	}
    }

    # moving window averaging
    foreach my $chr (@chr) {
	my @pos = (sort { $a <=> $b } keys %{$in{$chr}});
	my %average = ();
	INDEX:foreach my $index (0..$#pos) {
	    my $pos = $pos[$index];
	    my @vals = ($in{$chr}{$pos});
	    foreach my $extend (1..$window_size) {
		my $index2 = $index+$extend;
		last INDEX if ($index2 > $#pos);
		my $pos2 = $pos[$index2];
		last if ($pos2 > $pos+$window_size-1);
		push @vals, $in{$chr}{$pos2};
	    }
	    my $average = &mean(@vals);
	    my $mid_point = $pos + int($window_size / 2);
	    $average{$mid_point} = $average;
	}
	
	# detect chers
	@pos = (sort { $a <=> $b } keys %average);

	my $start = '';
	my $end = '';
	my $val = '';
	my @vals = ();
	foreach my $pos (@pos) {
	    $val = $average{$pos};
	    if ($val >= $min_score) {
		if ($start eq '') {
		    $start = $pos;
		    $end = $start;
		} else {
		    $end = $pos;
		}
		push @vals, $val;
	    } else {
		if ($start eq '') {
		} else {
		    my $average = &mean(@vals);
		    my $dist = $end - $start + 1;
		    if ($dist >= $min_width) {
			print "$chr\t$start\t$end\t$average\n";
		    }
		    @vals = ();
		    $start = '';
		    $end = '';
		}
	    }
	}
	unless ($end eq '') {
	    my $average = &mean(@vals);
	    my $dist = $end - $start + 1;
	    if ($dist >= $min_width) {
		print "$chr\t$start\t$end\t$average\n";
	    }
	}
    }

    close OUT;
}
		

sub mean {
    my @values = @_;
    my $sum = 0;
    foreach (@values) {
	$sum += $_;
    }
    return $sum / @values;
}
