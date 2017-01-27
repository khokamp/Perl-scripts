use warnings;
use strict;
use File::Basename;
use Getopt::Long;

my $usage = '

# script that converts mapped reads from BAM or SAM files 
# into strand-specific BigWig files
# normalisation can be achieved by specifying either 
# a reference size to which everything is scaled (-ref_size)
# or a file with scaling factors for each sample (-scale)
# requires UCSC\'s wigToBigWig to be installed and 
# a file with chromosome sizes (\'chrom.sizes\' by default)
# written by Karsten Hokamp (kahokamp@tcd.ie) in 2014/2015

';

my $chrom_sizes = 'chrom.sizes';
my $chrom_change = '';
my $neg = '';
my $limit = '';
my $no_flip = '';
my $keep_dups = '';
my $dups_only = '';
my $scaling_factor_file = '';
my $ref_size = '';
my $region = '';
my $do_not_overwrite = '';
my $help = '';
my $flip_strand = '';

&GetOptions(
    'help' => \$help,
    'flip_strand' => \$flip_strand,
    'do_not_overwrite' => \$do_not_overwrite,
    'region=s' => \$region,
    'ref_size=i' => \$ref_size,
    'dups_only' => \$dups_only,
    'keep_dups' => \$keep_dups,
    'chrom_change=s' => \$chrom_change,
    'chrom_sizes=s' => \$chrom_sizes,
    'neg' => \$neg,
    'no_flip' => \$no_flip,
    'limit=i' => \$limit,
    'scaling_factor_file|scale=s' => \$scaling_factor_file,
    );

if ($help) {
    print $usage;
    exit;
}

my %chrom_change = ();
if ($chrom_change) {
    my ($old, $new) = split /\t/, $chrom_change;
    $chrom_change{$old} = $new;
    print STDERR "Changing chromosome from '$old' to '$new'\n";
}


unless (-s $chrom_sizes) {
    die "No file '$chrom_sizes' found!\n";
} else {
    print STDERR "Using '$chrom_sizes' as file with chromosome sizes\n";
}

my %factor = ();
if ($scaling_factor_file) {
    open (IN, $scaling_factor_file)
        or die;
    while (<IN>) {
        chomp;
        my ($file, $factor, @rest) = split /\t/, $_;
        $factor{$file} = $factor;
    }
    close IN;
}

foreach my $file (@ARGV) {

    unless ($file =~ /(.*)\.bam/
	    or
	    $file =~ /(.*)\.sam/
	    or
	    $file =~ /(.*)\.bz2/) {
	warn "Not a bam, sam or .bz2 file: $file\n";
	next;
    }

    my $label = basename $1;

    my $time = time;

    print STDERR "Processing $label on ".(localtime)."...\n";

    my $factor = 1;
    if ($scaling_factor_file) {
        if (defined $factor{$file}) {
            $factor = 1/$factor{$file};
        } else {
            die "No factor for file '$file'";
	}
    } else {
	if ($ref_size) {
	    # get number of uniquely mapped reads in this file
	    my $size = 0;
	    print STDERR "Getting number of mapped reads: \n";
	    my $time = time;
	    if ($file =~ /\.bz2/) {	    
		open (IN, "bunzip2 -c \"$file\" |")
		    or die;
	    } else {
		if ($file =~ /\.sam$/) {
		    open (IN, "samtools view -S \"$file\" |")
			or die;
		} else {
		    open (IN, "samtools view \"$file\" |")
			or die;
		}
	    }
	    while (<IN>) {
		my @h = split /\t/, $_;
		my $flag = $h[1];
		
		# skip unmapped
		if ($flag & 4) {
		    next;
		}
		
		# skip reads that are paired but not mapped in proper pair:
		if ($flag & 1) {
		    unless ($flag & 2) {
			next;
		    }
		}
		
		# skip reads that are not primary alignment:
		if ($flag & 256) {
		    next;
		}
		
		# skip reads that are mapped to multiple locations
		if (defined $h[12]
		    and
		    $h[12] =~ /XS:i/) {
		    unless ($keep_dups
			    or
			    $dups_only) {
			next;
		    }
		}
		if ($dups_only) {
		    next unless (defined $h[12]
				 and
				 $h[12] =~ /XS:i/);
		}
		$size++;
	    }
	    close IN;

	    $factor = $ref_size / $size;
	    $time = time - $time;
	    print STDERR " $size ($time seconds)\n";
	}
    }
	    
    print STDERR "Adding each read as $factor\n";

    foreach my $sort_strand ('+', '-') {
#    foreach my $sort_strand ('+') {

	my $out = "$label.$sort_strand.wig";
	my $out2 = "$label.$sort_strand.bw";

	if (-s $out2) {
	    if ($do_not_overwrite) {
		warn "$out2 exists already - skipping\n";
		next;
	    }
	}

	my %skip = ();
	my %keep = ();

	if ($file =~ /\.bz2/) {	    
	    open (IN, "bunzip2 -c \"$file\" |")
		or die;
	} else {
	    if ($file =~ /\.sam$/) {
		open (IN, "samtools view -S \"$file\" $region |")
		    or die;
	    } else {
		open (IN, "samtools view \"$file\" $region |")
		    or die;
	    }
	}
	
	my $in = 0;
	my %map = ();
	while (<IN>) {
	    my @h = split /\t/, $_;
	    
	    $in++;

	    if ($limit ne '') {
		last if ($in > $limit);
	    }

	    if ($in % 1000000 == 0) {
		my $diff = time - $time;
		print STDERR " $in on ".(localtime)." ($diff)\n";
		$time = time;
	    }
	    
#	last if ($in > 10000);
	    
	    my $chr = $h[2];
	    
	    next if ($chr eq '*');

	    my $flag = $h[1];

	    # skip unmapped
	    if ($flag & 4) {
		$skip{$flag}++;
		next;
	    }

	    # skip reads that are paired but not mapped in proper pair:
	    if ($flag & 1) {
		unless ($flag & 2) {
		    $skip{$flag}++;
		    next;
		}
	    }

	    # skip reads that are not primary alignment:
	    if ($flag & 256) {
		$skip{$flag}++;
		next;
	    }

	    # skip reads that are mapped to multiple locations
	    if (defined $h[12]
		and
		$h[12] =~ /XS:i/) {
		unless ($keep_dups
			or
			$dups_only) {
		    $skip{'multip-map'}++;
		    next;
		}
	    }
	    if ($dups_only) {
		next unless (defined $h[12]
			     and
			     $h[12] =~ /XS:i/);
	    }
	    
	    
	    my $strand = $flag;
	    if ($strand & 16 ) {
		$strand = '-';
		if ($flip_strand) {
		    $strand = '+';
		}
	    } else {
		$strand = '+';
		if ($flip_strand) {
		    $strand = '-';
		}
	    }

	    if ($flag & 128) {
		# second in pair
		# reverse strand
		unless ($no_flip) {
		    if ($strand eq '-') {
			$strand = '+';
			if ($flip_strand) {
			    $strand = '-';
			}
		    } else {
			$strand = '-';
			if ($flip_strand) {
			    $strand = '+';
			}
		    }
		}
	    }

	    unless ($strand eq $sort_strand) {
		$skip{$flag}++;
		next;
	    }
	    
	    $keep{$flag}++;

	    my $start = $h[3];
	    my $temp = length($h[9]);
	    my $end = $start + $temp - 1;

	    my $cigar = $h[5];
	    if ($cigar =~ /^(\d+)S/) {
#		$start += $1;
		# start position is already corrected (in Bowtie2)
		# but length of match must be shortened
		$end -= $1;
	    }
	    if ($cigar =~ /(\d+)S$/) {
		$end -= $1;
	    }
	    
	    foreach ($start..$end) {
		$map{$chr}{$_} += $factor;
	    }
	}
	close IN;
	
	print STDERR "\nFinished reading on ".(localtime)."...\n";
	
	print STDERR "\nKept the following flags:\n";
	foreach my $flag (sort { $a <=> $b } keys %keep) {
	    print STDERR "$flag ($keep{$flag})\n";
	}
	print STDERR "\nSkipped the following flags:\n";
	foreach my $flag (sort { $skip{$a} <=> $skip{$b} } keys %skip) {
	    print STDERR "$flag ($skip{$flag})\n";
	}

	open (OUT, ">$out")
	    or die;
	
	foreach my $chr (sort keys %map) {
	    
	    my $diff = time - $time;
	    print STDERR " $chr ($sort_strand) on ".(localtime)." ($diff)\n";
	    $time = time;
	    
	    # write wig file:
	    my $chr_out = $chr;
	    if (defined $chrom_change{$chr}) {
		$chr_out = $chrom_change{$chr};
	    }
	    print OUT "variableStep chrom=$chr_out span=1\n";
	    
	    foreach my $pos (sort { $a <=> $b } keys %{$map{$chr}}) {
		my $val = $map{$chr}{$pos};
		if ($neg
		    and
		    $sort_strand eq '-') {
		    $val = -1 * $val;
		}
		print OUT "$pos\t$val\n";
	    }
	}
	
	close OUT;
	
	my $job = "wigToBigWig -clip \"$out\" $chrom_sizes \"$out2\" 2> /tmp/err.$$";
	print STDERR "job: $job\n";
	system "$job";
	
	print STDERR "Changed to BigWig $out2 on ".(localtime)."...\n";
	
    }
}

