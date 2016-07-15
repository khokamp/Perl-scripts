use warnings;
use strict;
use Getopt::Long;

# convert a GenBank file into GFF format
# run like this:
# perl gb2gff3.pl file.gb
# output is written to file.gff3
# use -add_annotation flag to include gene name and product description in the output

my $add_annotation = '';

&GetOptions(
    'add_annotation' => \$add_annotation,
    );

foreach my $file (@ARGV) {
    my $name = $file;    
    $name =~ s/(.*?)\..*/$1/;
    warn "Processing $file...\n";

    my %annotation = ();
    if ($add_annotation) {
	# get annotation first
	open (IN, $file)
	    or die;
	my $type = '';
	my $cds = '';
	while (<IN>) {
	    last if (/^FEATURES/);
	}

	while (<IN>) {
	    chomp;
	    if (/^\s{4,}(\w+)\s{4,}(\w.+)/) {
		$type = $1;
		$cds = $2;
	    } else {
		if (/^\s+\/(.+)="(.+)/) {
		    my $cat = $1;
		    next unless ($cat eq 'gene'
				or
				$cat eq 'product');
		    my $val = $2;
		    until ($val =~ /"$/) {
			$_ = <IN>;
			chomp;
			s/^\s+/ /;
			$val .= $_;
		    }
		    $val =~ s/\"$//;
		    $annotation{$cds}{$type}{$cat} = $val;
		}
	    }
	}
	close IN;
    }
	    
    open (IN, $file)
	or die;
    open (OUT, ">$name.gff3")
	or die;

    my @out = ();
    my %pos = ();
    my $num = 0;
    my $cds = '';
    while (<IN>) {
	last if (/^FEATURES/);
    }

    while (<IN>) {
	last if (/^ORIGIN/);
	chomp;
	s/\cM//;
	if (/systematic_id="(.*)"/
	    or
	    /locus_tag="(.*)"/) {
	    my $id = $1;
	    $num++;
	    my $part = '';
	    my $add = '';
	    if (@out > 1) {
		$add = ' part ';
		$part = 1;
	    }
	    my $ann = '';
	    if ($add_annotation) {
		my @ann = ();
		if (defined $annotation{$cds}{'CDS'}{'gene'}) {
		    push @ann, "Name=$annotation{$cds}{'CDS'}{'gene'}";
		}
		if (defined $annotation{$cds}{'CDS'}{'product'}) {
                    push @ann, "Description=$annotation{$cds}{'CDS'}{'product'}";
                }
		$ann = ';' . (join ';', @ann);
	    }

	    foreach (@out) {
		my ($start, $end, $strand) = @$_;
		if (defined $pos{$start}{$end}) {
		    die "$start..$end defined more than once!\n";
		}
		$pos{$start}{$end} = "$name\tpredicted\tORF\t$start\t$end\t.\t$strand\t.\tID=$id$add$part$ann\n";
		$part++;
	    }
	    @out = ();
	} elsif (/   CDS\s+(.*)/) {
	    $cds = $1;
	    my $strand = '+';
	    my ($start, $end);
	    if ($cds =~ /^join\((.*)\)$/
		or
		$cds =~ /^complement\(join\((.*)\)\)$/
		) {
		my $cds2 = $1;
		if ($cds =~ /complement/) {
		    $strand = '-';
		}
		foreach (split /,/, $cds2) { #/) {
		    unless (/^(\d+)\.\.(\d+)$/) {
			die "unknown part in join: '$_'\n";
		    }
		    $start = $1;
		    $end = $2;
		    push @out, [ $start, $end, $strand ];
		}
	    } elsif ($cds =~ /^complement\((\d+)\.\.(\d+)\)$/) {
		$start = $1;
		$end = $2;
		$strand = '-';
		@out = [$start, $end, $strand];
	    } elsif ($cds =~ /^(\d+)\.\.(\d+)$/) {
		$start = $1;
		$end = $2;
		$strand = '+';
		@out = [$start, $end, $strand];
	    } else {
		die "Unknown format: '$cds'!\n";
	    }
	}
    }
    
    close IN;

    warn "$num entries extracted\n";

#    if (/SQ\s+Sequence (\d+)/) {
#	print OUT "##gff-version 3
###Note: See http://song.sourceforge.net
#
#$name	predicted	contig	1	$1	.	.	.	Name=$name
#";
#
#    } else {
#	warn "Unknown first line: $_\n";
#    }
    foreach my $start (sort { $a <=> $b } keys %pos) {
	foreach my $end (sort { $a <=> $b } keys %{$pos{$start}}) {
	    print OUT $pos{$start}{$end};
	}
    }
}
    
