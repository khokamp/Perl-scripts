#!usr/bin/perl

# SCRIPT DETAILS:-
# Perl script to calculate the number of reads mapped to each annotated genes.
# Syntax: perl reads_count_gff.pl annotation.gff mapped_reads1.bam [mapped_reads2.bam ...] > read_counts.txt

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my $annotation_file = '';
my $add_position = '';
my $out_dir = '.';
my $add_antisense = '';
my $no_flip = '';
my $keep_dups = '';
my $dups_only = '';
my $annotation_style = '';

&GetOptions(
    'dups_only' => \$dups_only,
    'keep_dups' => \$keep_dups,
    'add_antisense' => \$add_antisense,
    'out_dir=s' => \$out_dir,
    'annotation_file=s' => \$annotation_file,
    'annotation_style=s' => \$annotation_style,
    'add_position=i' => \$add_position,
    );

open (IN, $annotation_file)
    or die;
my %overlap = ();
my %ann = ();
my %rna = ();

# Annotation file might look like this:
#ST4-74.fa       Kroeger_2013    CDS     169     255     .       +       .       ID=SL1344_0001;Name=thrL;description=thr operon leader peptide

while (<IN>) {
    next unless (/\w/);
    next if (/^#/);
    chomp;
    s/\cM//;
    my $id = '';
    my ($chr, $source, $type, $start, $end, $undef1, $strand, $undef2, $info) = split /\t/, $_;
    next if ($type eq 'contig');
    if ($annotation_style eq 'tab') {
	($id, $chr, $strand, $start, $end) = split /\t/, $_;
	$info = "ID=$id;";
    }
    if ($info =~ /ID=(.+?)\;/) {
	$id = $1;
    } else {
	# skip dummy entry
	next if ($info =~ /XXX/);
	die "No id: $info\n";
    }

    my $name = '';
    if ($info =~ /Name=(.+?);/) {
	$name = $1;
    }

    if ($start > $end) {
	($start, $end) = ($end, $start);
    }

    if ($add_position ne '') {
	if ($start >= $add_position) {
	    $start++;
	    $end++;
	}
    }

    if (defined $name) {
	$ann{$id}{'name'} = $name if ($name =~ /\w/);
	if ($name =~ /\w/) {
	    $ann{$id}{'name'} = $name;
	} else {
	    $ann{$id}{'name'} = $id;
	}
    }

    $ann{$id}{'strand'} = $strand;
    $ann{$id}{'start'} = $start;
    $ann{$id}{'end'} = $end;
    $ann{$id}{'chr'} = $chr;

#    my $id2 = "$id";
#    if ($name = /\w/) {
#	$id2 .= ", $name";
#    }
    my $id2 = $id;

    foreach ($start..$end) {
	if (defined $rna{$chr}{$strand}{$_}) {
	    $overlap{$id2}{$rna{$chr}{$strand}{$_}}++;
	}
	$rna{$chr}{$strand}{$_} = $id;
    }

    if ($add_antisense) {
	if ($strand eq '+') {
	    $strand = '-';
	} else {
	    $strand = '+';
	}

	$id .= '_as';
	$id2 .= '_as';

	$ann{$id}{'strand'} = $strand;
	$ann{$id}{'start'} = $start;
	$ann{$id}{'end'} = $end;
	$ann{$id}{'chr'} = $chr;

	my $id2 = $id;

	foreach ($start..$end) {
	    if (defined $rna{$chr}{$strand}{$_}) {
		$overlap{$id2}{$rna{$chr}{$strand}{$_}}++;
	    }
	    $rna{$chr}{$strand}{$_} = $id;
	}
    }
}
close IN;

# report overlap amongst annotated features:
print "Overlapping annotation:\n";
my %overlap2 = ();
foreach my $id (sort keys %overlap) {
    my @ids = sort keys %{$overlap{$id}};
    foreach my $id2 (@ids) {
	my $len = $overlap{$id}{$id2};
	$overlap2{$len}{"$id overlaps $id2"}++;
    }
}

foreach my $len (sort { $a <=> $b } keys %overlap2) {
    foreach my $ids (sort keys %{$overlap2{$len}}) {
	print "Overlap: $ids by $len bp\n";
    }
}

# process input files
foreach my $file (@ARGV) {

    print STDERR "Processing $file on ".(localtime)." ... ";
    my $time = time;

    my ($tag, $suffix);
    my $base = basename($file);
    if ($base =~ /^(\w+?)\.(.+)/) {
	$tag = $1;
	$suffix = $2;
    } elsif ($base =~ /^(.+)\.sorted\.bam/) {
	$tag = $1;
    } else {
	die "Unknown format: '$file'!";
    }

    my $out = "$out_dir/$tag.count";
    
    if (-e $out) {
	die "Outfile '$out' exists already!";
    }

    open (OUT, ">$out")
	or die "Can't write to $out: $!";

    if ($file =~ /\.bam$/) {
	open (READ, "samtools view \"$file\" |")
	    or die;
    } elsif ($file =~ /\.bz2/) {
	open (READ, "bunzip2 -c \"$file\" |")
	    or die;
    } else {
	die "Invalid file ending: $file\n";
    }

    my %cover = ();
    my %stats = ();
    my %skip = ();
    my %reads = ();
    while(<READ>){
	$stats{'reads'}++;
	chomp;
	my @element = split /\t/, $_, 13;
	my $strand = $element[1];
	my $start = $element[3];
	$reads{$start}++;
	my $temp = length($element[9]);
	my $end = $start + $temp - 1;
	my $chr = $element[2];
	my $cigar = $element[5];
	if ($cigar =~ /^(\d+)S/) {
#	    $start += $1;
	    # start position is already corrected (in Bowtie2)
	    # but length of match must be shortened
	    $end -= $1;
	}
	if ($cigar =~ /(\d+)S$/) {
	    $end -= $1;
	}
	my $len = $end - $start + 1;
	$stats{'coverage'} += $len;
	$stats{"lengths $len"}++;
#	if ($strand eq '0') {
#	    $strand = '+';
#	} else {
#	    $strand = '-';
#	}

	my $flag = $strand;
	
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
	if (defined $element[12]
	    and
	    $element[12] =~ /XS:i/) {
	    unless ($keep_dups
		    or
		    $dups_only) {
		$skip{'multip-map'}++;
		next;
	    }
	}
	if ($dups_only) {
	    next unless (defined $element[12]
			 and
			 $element[12] =~ /XS:i/);
	}

	# assign strand
	if ($flag & 16 ) {
	    $strand = '-';
	} else {
	    $strand = '+';
	}
	
	if ($flag & 128) {
	    # second in pair
	    # reverse strand
	    unless ($no_flip) {
		if ($strand eq '-') {
		    $strand = '+';
		} else {
		    $strand = '-';
		}
	    }
	}
	
	
	$stats{"$strand strand"}++;

	$stats{"chr $chr"}++;

	my $start_id = '';
	my $end_id = '';
	if (defined $rna{$chr}{$strand}{$start}) {
	    $start_id = $rna{$chr}{$strand}{$start};
	    $stats{'cover start'}++;
	}

	if (defined $rna{$chr}{$strand}{$end}) {
	    $end_id = $rna{$chr}{$strand}{$end};
	    $stats{'cover end'}++;
	}

	my $id = '';
	if ($start_id and $end_id) {
	    $stats{'cover'}++;
	    if ($start_id ne $end_id) {
		if ($strand eq '+') {
		    $id = $start_id;
		    $stats{'covering different features - assigning start'}++;
		} else {
		    $id = $end_id;
		    $stats{'covering different features - assigning end'}++;
		}
	    } else {
		$id = $start_id;
	    }
	} else {
	    if ($start_id) {
		$stats{'cover'}++;
		$id = $start_id;
	    } elsif ($end_id) {
		$stats{'cover'}++;
		$id = $end_id;
	    } else {
		$stats{'no cover'}++;
	    }
	}

	if ($id) {
	    $cover{$id}++;
	}
    }
    close READ;

    print STDERR "\nSkipped the following flags:\n";
    foreach my $flag (sort { $skip{$a} <=> $skip{$b} } keys %skip) {
	print STDERR "$flag ($skip{$flag})\n";
    }

    my $duration = time - $time;
    print STDERR " done in $duration seconds\n";

    foreach my $id (sort keys %ann){
	my $val = 0;
	if (defined $cover{$id}) {
	    $val = $cover{$id};
	}
	print OUT "$id\t$val\n";
	while (length($val) < 5) {
	    $val = ' '.$val;
	}
	$stats{"coverage: $val"}++;
    }
    close OUT;

    $stats{'average read length'} = $stats{'coverage'} / $stats{'reads'};

    print "\nStats for $file:\n";
    foreach my $key (sort keys %stats) {
	print "Stats: $key -> $stats{$key}\n";
    }

    print "\n";
}

