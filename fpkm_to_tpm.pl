use warnings;
use strict;

# Karsten Hokamp, Trinity College Dublin, 2019
# turn RPKM or FPKM values into TPM values
# this currently works for files that have gene names in the first column
# for each subsequent column the values are summed up and used in the formula 
# taken from here: https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
# to calculate TPM values

my $header = <>;
print $header;

my %in = ();
my @in = ();

# go through all lines to collect values and sum up columns
while (<>) {
    chomp;
    push @in, $_;
    my @h = split /\t/, $_;
    foreach (1..$#h) {
	$in{$_} += $h[$_];
    }
}

# go through collected values to divide by column sum and multiple by a million:
foreach (@in) {
    my @h = split /\t/, $_;
    foreach (1..$#h) {
	$h[$_] = $h[$_] / $in{$_} * 1000000;
    }
    print "".(join "\t", @h)."\n";
}

    


