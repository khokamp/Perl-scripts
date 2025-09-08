use warnings;
use strict;
use Getopt::Long;

# Karsten Hokamp, Trinity College Dublin, 2024
# filter gene count matrix by TPM
# by default only keeps genes that are expressed with at least 2 TPM in more than half of the samples within at least one of the conditions
# (conditions are recognised by start of column headers ending in TPM)

my $min_tpm = 2;
my $exclude = '';
my $rule = 'more than half';

&GetOptions(
    'exclude=s' => \$exclude,
    'min_tpm=s' => \$min_tpm,
    'rule=s' => \$rule,
    );

my $header = <>;
$header =~ s/\.bam//g;
print $header;
chomp $header;
my @header = split /\t/, $header;
my %header = ();
my %groups = ();

my $i = 0;
foreach my $head (@header) {
    $header{$head} = $i;
    $i++;
    if ($head =~ /(.+)_\d+ TPM/
	or
	$head =~ /(.+)\d+ TPM/) {
	my $group = $1;
	if ($group eq $exclude) {
	    warn "Skipping $head\n";
	    next;
	}
	$groups{$group}{$head}++;
    }
}

foreach my $group (sort keys %groups) {
    my @files = sort keys %{$groups{$group}};
    my $num = scalar @files;
    print STDERR "Group $group - $num files: @files\n";
}

my $kept = 0;
my $skipped = 0;
while (<>) {
    chomp;
    my @h = split /\t/, $_;
    my $keep = 0;
    foreach my $group (sort keys %groups) {
	my $max = 0;
	my @files = sort keys %{$groups{$group}};
	my $num = scalar @files;
	my $threshold = int($num / 2);
	foreach my $sample (@files) {
	    my $col = $header{$sample};
	    my $val = $h[$col];
	    if ($val >= $min_tpm) {
		$max++;
	    }
	}
	if ($rule eq 'more than half') {
	    if ($max > $threshold) {
		$keep = 1;
	    }
	} elsif ($rule eq 'at least half') {
	    if ($max >= $threshold) {
                $keep = 1;
            }
	} else {
	    die "Unknown rule: '$rule'";
	}
	last if ($keep);
    }
    if ($keep) {
	print "$_\n";
	$kept++;
    } else {
	$skipped++;
    }
}

print STDERR "Printed $kept entries, skipped $skipped (no group had more than half of samples with TPM >= $min_tpm)\n";
		
