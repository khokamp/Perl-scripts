use warnings;
use strict;

# Karsten Hokamp, Trinity College Dublin, 2015
# Perl script that generates BigWig files from RMAT output
# requires file with chromosome sizes and wigToBigWig tool from UCSC
# input consists of BED files containing chr, start and end position of a probe,
# as well as the associated MAT score (see example lines below)
# output will be written to files with the same name as input files but .bw instead of .bed ending                                                                                                      
# values from overlapping probes are averaged using the mean

my $chrom_sizes = 'chrom_sizes.txt';

unless (-s $chrom_sizes) {
    die "No file '$chrom_sizes' found!";
}

my %chr = ();
my $chr = 0;
foreach ('I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI') {
    $chr++;
    $chr{$chr} = 'chr'.$_;
}

foreach my $file (@ARGV) {

    print "Processing '$file'...\n";
    
    unless ($file =~ /\.bed/) {
	warn "Invalid file name '$file'\n";
	next;
    }

    open (IN, $file)
	or die;
    my $out = $file;    
    $out =~ s/\.bed/\.wig/;

    my $out2 = $file;    
    $out2 =~ s/\.bed/\.bw/;

    open (OUT, ">$out")
	or die;
    
#chrI    0       24      -4.70441136395268
#chrI    4       28      -4.82038226036296
#chrI    8       32      -4.97036782313776

    my %cover = ();
    while (<IN>) {
	next unless (/\w/);
	chomp;
	s/\cM//;
	my ($chr, $start, $end, $score) = split /\t/, $_;
	$start++;
	$end++;
	if ($chr =~ /^\"chr(\d+)\"/) {
	    $chr = $1;
	    if (defined $chr{$chr}) {
		$chr = $chr{$chr};
	    }
	} elsif ($chr =~ /(chr\w+)/) {
	    $chr = $1;
	} else {
	    die "Unknown chr: '$chr'\n";
	}

	foreach ($start .. $end) {
	    push @{$cover{$chr}{$_}}, $score;
	}
    }

    foreach my $chr (sort keys %cover) {
	print OUT "variableStep chrom=$chr span=1\n";
	foreach my $pos (sort { $a <=> $b } keys %{$cover{$chr}}) {
	    my $val = &mean(@{$cover{$chr}{$pos}});	    	
	    print OUT "$pos\t$val\n";
	}
    }
    close OUT;
    close IN;

    system "wigToBigWig -clip \"$out\" \"$chrom_sizes\" \"$out2\"";
   
    print STDERR "Changed to BigWig $out2 on ".(localtime)."...\n";
}

	

sub mean {
    my @values = @_;
    my $sum = 0;
    foreach (@values) {
	$sum += $_;
    }
    return $sum / @values;
}


