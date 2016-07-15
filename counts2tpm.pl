use warnings;
use strict;
use Getopt::Long;

# Karsten Hokamp, Trinity College Dublin, 2013
# calculate various metrics from read counts summarised from RNA-seq data
# TPM is the implementation of Li et al (Bioinformatics, 2010, doi:10.1093/bioinformatics/btp692)
# TPM_Wagner is the implementation of Wagner et al (Theory Biosci., 2012, doi:10.1007/s12064-012-0162-3)
# RPKM as defined by Mortazavi et al (Nature Methods, 2008)
# 
# Note: TPM and TPM_Wagner are different ways of getting to the same results!

my $USAGE = "
$0 [OPTIONS] data_file > output

data_file contains first column with geneID followed by 
one or more columns with read counts for each library

OPTIONS:
-method           method to apply (default 'TPM')
-annotation_file  file with geneID and annootation info
-gene_length_file file with geneID and transcript length
-rounding_ratio   number of digits after decimal point for ratios
-rounding_abs     number of digits after decimal point for non-ratios
-read_length      average read length (only used for generation of 
                  scaling factors in additional output)
-help             this help text

";

my $method = 'TPM';
my $annotation_file = '';
my $gene_length_file = '';
my $rounding_ratio = '';
my $rounding_abs = '';
my $read_length = 50;
my $help = '';
my $skip = '';
my $length_in_data_file = '';
my $replicates = '';

&GetOptions(
    'read_length=i' => \$read_length,
    'method=s' => \$method,
    'skip=s' => \$skip,
    'rounding_ratio|rounding=s' => \$rounding_ratio,
    'rounding_abs=s' => \$rounding_abs,
    'annotation_file=s' => \$annotation_file,
    'gene_length_file|lengths=s' => \$gene_length_file,
    'length_in_data_file' => \$length_in_data_file,
    'replicates=i' => \$replicates,
    'help' => \$help,
    );

if ($help) {
    print $USAGE;
    exit;
}

my %skip_lib_num = ();
foreach (split /,/, $skip) { #/) {
    $skip_lib_num{$_}++;
}

my @methods = ('raw', 'RPKM', 'TPM', 'TPM_Wagner');
if ($method eq 'all') {
    $method = join ',', (@methods);
}

if ($method =~ /P/) {
    unless ($gene_length_file) {
	die "Gene length file needed!\n" unless ($length_in_data_file);
    }
}

# read in annotation (optional)
my %ann = ();
my @ann = ();
my $ann_size = 0;
if ($annotation_file) {
    open (IN, $annotation_file)
	or die;
    my $header = <IN>;
    $header =~ s/\cM//;
    chomp $header;
    @ann = split /\t/, $header;
    shift @ann;
    while (<IN>) {
	chomp;
	s/\cM//;
	my @h = split /\t/, $_;
	my $id = shift @h;
	$ann{$id} = join "\t", @h;
	if (@h > $ann_size) {
	    $ann_size = @h;
	}
    }
    close IN;
}

# read in gene lengths (required for TPM and RPKM)
my %lengths = ();
if ($gene_length_file or $length_in_data_file) {
    if ($length_in_data_file) {
	open (IN, $ARGV[0])
	    or die;
    } else {
	open (IN, $gene_length_file)
	    or die "Can't open gene length file '$gene_length_file': $!";
    }
    while (<IN>) {
	chomp;
	s/\cM//;
	my ($gene, $len, @rest) = split /\t/, $_;
	$lengths{$gene} = $len;
    }
    close IN;
}

# read in main data file with read counts per gene:
my $header = <>;
chomp $header;
$header =~ s/\cM//;
my %headers = &get_headers($header);

my @all_libs = split /\t/, $header;
# first column contains gene IDs
shift @all_libs;
if ($length_in_data_file) {
    shift @all_libs;
}


# read in raw counts
my %in = ();
my %raw_lib_size = ();
my %skip_lib = ();
while (<>) {
    chomp;
    s/\cM//;
    my @h = split /\t/, $_;
    my $id = $h[0];
    $id =~ s/\"//g;

    foreach my $lib (sort @all_libs) {
	unless (defined $headers{$lib}) {
	    die;
	}
	my $col = $headers{$lib};
	if (defined $skip_lib_num{$col}) {
	    $skip_lib{$lib}++;
	    next;
	}
	my $val = $h[$col];
	$in{'raw'}{$id}{$lib} = $val;
	$raw_lib_size{$lib} += $val;
    }
}

# calculate factor for TPM_Wagner if required
my %skip = ();
my %T_W = ();
if ($method =~ /TPM_Wagner/) {
    foreach my $id (sort keys %{$in{'raw'}}) {
	unless (defined $lengths{$id}) {
	    warn "No gene length for $id!";
	    $skip{$id}++;
	    next;
	}
	my $len = $lengths{$id};
	
	foreach my $lib (@all_libs) {
	    next if (defined $skip_lib{$lib});
	    my $val = $in{'raw'}{$id}{$lib};
	    $T_W{$lib} += $val / $len;
	}
    }
}

# calculate factor for TPM if required
my %T_L = ();
if ($method =~ /TPM/) {
    foreach my $id (sort keys %{$in{'raw'}}) {
	unless (defined $lengths{$id}) {
	    warn "No gene length for $id!\n";
	    $skip{$id}++;
	    next;
	}
	my $len = $lengths{$id};
	if ($len == 0) {
	    die "Length is zero for $id\n";
	}
	foreach my $lib (@all_libs) {
	    next if (defined $skip_lib{$lib});
	    my $val = $in{'raw'}{$id}{$lib} / $raw_lib_size{$lib};
	    $T_L{$lib} += $val / $len;
	}
    }
}


# carry out selected methods of normalisation
my %out = ();
foreach my $type (split /,/, $method) { #/) {
    &normalise($type);
}


# calculate TPM size factors for normalisation of raw libraries:
# (could be used to adjust bigwig tracks in Jbrowse or gr files for IGB)
if (defined $in{'TPM'}) {
    my $tpm_factor_file = 'tpm_factors.txt';
    open (TPM, ">$tpm_factor_file")
	or die;
    print STDERR "Printing TPM factors to $tpm_factor_file\n";
    foreach my $lib (sort keys %raw_lib_size) {
	my $raw = $raw_lib_size{$lib};
	my $tpm = 0;
	foreach my $gene (keys %{$in{'TPM'}}) {
	    my $val = $in{'TPM'}{$gene}{$lib};
	    my $len = $lengths{$gene};
	    $tpm += $val * $len;
	}
	$tpm /= $read_length;
	my $factor = $tpm / $raw;
	print TPM "$lib\t$factor\n";
    }
    close TPM;
}


# print output:
# 1. gene, annotation and length
# 2. all raw and derived values (rpkm, tpm)
# 3. the ratios for each method
my @header = ('gene', @ann);
if ($gene_length_file
    or
    $length_in_data_file) {
    push @header, 'gene length';
}
foreach (@methods) {
    my %done = ();
    if (defined $in{$_}) {
	my $counter = 0; 
	foreach my $lib (@all_libs) {
	    next if (defined $done{$lib});
	    next if (defined $skip_lib{$lib});
	    push @header, "$lib $_";
	    $done{$lib}++;
	    $counter++;
	    if ($replicates) {
		if ($counter % $replicates == 0) {
		    push @header, "merge $_";
		}
	    }		
	}
	push @header, "min $_\tmax $_\tmean $_";
    }
}

print "".(join "\t", @header)."\n";

foreach my $gene (sort keys %{$in{'raw'}}) {
    next if (defined $skip{$gene});
    my @out = ($gene);
    if (%ann) {
	if (defined $ann{$gene}) {
	    push @out, $ann{$gene};
	} else {
	    my $ann = join "\t", ('' x $ann_size);
	    push @out, $ann;
	}	    
    }
    if (%lengths) {
	push @out, $lengths{$gene};
    }

    foreach my $type (@methods) {
	my %done = ();
	if (defined $in{$type}) {
	    my @vals = ();
	    foreach my $lib (@all_libs) {
		next if (defined $done{$lib});
		next if (defined $skip_lib{$lib});
		push @out, $in{$type}{$gene}{$lib} unless ($replicates);
		push @vals, $in{$type}{$gene}{$lib};
		$done{$lib}++;
	    }
	    if ($replicates) {
		# merge replicates
		my @vals2 = ();
		while (@vals) {
		    my @reps = ();
		    foreach (1..$replicates) {
			my $val = shift @vals;
			push @reps, $val;
		    }
		    my $mean = &mean(@reps);
		    push @vals2, $mean;
		    push @out, @reps, $mean;
		}
		@vals = @vals2;
	    }
	    @vals = sort { $a <=> $b } @vals;
	    my $mean = sprintf "%.1f", &mean(@vals);
	    push @out, "$vals[0]\t$vals[-1]\t$mean";
	}
    }

    
    unless (defined $out[1]) {
	die;
    }
    print "".(join "\t", @out)."\n";
}


################
# Subroutines: #
################

	    
sub normalise {
    my $method = shift;
    
    if ($method eq 'RPKM') {
	foreach my $gene (sort keys %{$in{'raw'}}) {
	    next if (defined $skip{$gene});
	    unless (defined $lengths{$gene}) {
		warn "No length for $gene\n";
		$skip{$gene}++;
		next;
	    }
	    my $len = $lengths{$gene};
	    foreach my $lib (@all_libs) {
		next if (defined $skip_lib{$lib});
		my $val = $in{'raw'}{$gene}{$lib};
		my $derived = $val / ($len/1000) / ($raw_lib_size{$lib}/1000000);
		if ($rounding_abs) {
		    $derived = sprintf "%.${rounding_abs}f", $derived;
		    $derived =~ s/\.0+$//;
		}
		$in{$method}{$gene}{$lib} = $derived;
	    }
	}

    } elsif ($method eq 'TPM_Wagner') {
	foreach my $gene (sort keys %{$in{'raw'}}) {
	    next if (defined $skip{$gene});
	    unless (defined $lengths{$gene}) {
		warn "No length for $gene\n";
                $skip{$gene}++;
                next;
	    }
	    my $len = $lengths{$gene};
	    foreach my $lib (@all_libs) {
		next if (defined $skip_lib{$lib});
		my $val = $in{'raw'}{$gene}{$lib};
		my $derived = $val * 1000000 / ($len * $T_W{$lib});
		if ($rounding_abs) {
		    $derived = sprintf "%.${rounding_abs}f", $derived;
		    $derived =~ s/\.0+$//;
		}
		$in{$method}{$gene}{$lib} = $derived;
	    }
	}

    } elsif ($method eq 'TPM') {
	foreach my $gene (sort keys %{$in{'raw'}}) {
	    next if (defined $skip{$gene});
	    unless (defined $lengths{$gene}) {
		warn "No length for $gene\n";
                $skip{$gene}++;
                next;
	    }
	    my $len = $lengths{$gene};
	    foreach my $lib (@all_libs) {
		next if (defined $skip_lib{$lib});
		my $val = $in{'raw'}{$gene}{$lib};
		my $derived = ($val / $raw_lib_size{$lib}) * 1000000 / ($len * $T_L{$lib});
		if ($rounding_abs) {
		    $derived = sprintf "%.${rounding_abs}f", $derived;
		    $derived =~ s/\.0+$//;
		}
		$in{$method}{$gene}{$lib} = $derived;
	    }
	}
    }
}


sub mean {
    my @values = @_;
    my $sum = 0;
    foreach (@values) {
	$sum += $_;
    }
    return $sum / @values;
}


sub get_headers {
    my %out = ();
    my $header = shift;
    chomp $header;
    my $count = 0;
    my @tmp = split /\t/, $header;
    foreach (@tmp) {
	$out{$_} = $count;
	$count++;
    }
    return %out;
}
