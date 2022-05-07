#!/usr/bin/perl
use strict;
use warnings;
use 5.010;
use Spreadsheet::Read qw(ReadData);
use Getopt::Long;
use Excel::Writer::XLSX;

# read data from Easychair (one large Excel file)
# and combine data from several sheets into one table
# which is written as an Excel file
# (originally implemented for BCC2020)
    
#    Copyright (C) 2020  Karsten Hokamp
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see https://www.gnu.org/licenses/.

my $submissions_file = '';
my $out_file = '';
my $overwrite = '';

&GetOptions(
    'submissions_file=s' => \$submissions_file,
    'output|out_file=s' => \$out_file,
    'overwrite' => \$overwrite,
    );

if ($submissions_file eq '') {
    die "Please provide the name of a data file downloaded from Easychair via the '-submissions' parameter\n";
}

print "Processing $submissions_file\n";

if ($out_file eq '') {
    # name output after submission file:
    $out_file = $submissions_file;
    $out_file =~ s/\..+//;
    $out_file .= '_mashup.xlsx';
}
print "Writing output to $out_file\n";

if (-e $out_file) {
    warn "WARNING: $out_file exists already\n";
    unless ($overwrite) {
	die "Need -overwrite option to overwrite the file\n";
    } else {
	warn "file will be overwritten\n";
    }
}

# Create a new Excel workbook
my $workbook  = Excel::Writer::XLSX->new( $out_file );
my $worksheet = $workbook->add_worksheet();
my $formatw = $workbook->add_format();
$formatw->set_text_wrap();

my $formatb = $workbook->add_format();
$formatb->set_bold();

my $formatr = $workbook->add_format(align  => 'right');
#my $formatc = $workbook->add_format(align  => 'center');
my $formatc = $workbook->add_format(align  => 'top');


# read submission id and title from Easychair (copied into Excel):
#my $book = ReadData($submissions_file);
my $book = Spreadsheet::Read->new ($submissions_file);

my @sheets_found = $book->sheets;
my %sheets_found = ();
foreach (@sheets_found) {
    $sheets_found{$_}++;
}
# need the following sheets:
# Submissions: id and title
# Submission field values: requested presentation format
# Reviews: total score
# Review texts: reviewers' comments
# Review scores: individual scores
# Review form fields: names for field numbers

my %input_fields = (
    'Submissions', [ '#', 'title', 'authors', 'decision'],
    'Reviews', [ 'submission #', '#', 'total score' ],
    'Review texts', [ 'review #', 'field #', 'text' ],
    'Review scores', [ 'review #', 'field #', 'score' ],
    'Submission field values', [  'submission #', 'field name', 'value' ],
    'Review form fields', [ '#', 'title', 'has scores?', 'has text?' ],
    );

my @sheets = (sort keys %input_fields);

foreach my $sheet (@sheets) {
    unless (defined $sheets_found{$sheet}) {
	die "No sheet '$sheet' in $submissions_file\n";
    }
}


# read info
my %input = ();
foreach my $sheet (@sheets) {
    my $sheet_in = $book->sheet ($sheet);
    my @rows = $sheet_in->rows ();

    my %headers = ();
    foreach my $j (0 .. $#{$rows[0]}) {
	my $val = $rows[0][$j];
	$headers{$val} = $j;
    }

    my @headers = @{$input_fields{$sheet}};
    foreach my $header (@headers) {
	unless (defined $headers{$header}) {
	    die "no column with header '$header' found in sheet $sheet\n";
	}
    }

    # extract data
    my $ref = shift @headers;
    foreach my $i (1 .. $#rows) {
	# get the reference
	my $col = $headers{$ref};
	my $ref_val = $rows[$i][$col];

	# associate other values with it
	my @vals = ();
	foreach my $header (@headers) {
	    my $col = $headers{$header};
	    my $val = $rows[$i][$col];
	    push @vals, $val;
	}

	# store records of varying complexity (not the most elegant solution, but works):
	if (@headers == 1) {
	    $input{$sheet}{$ref_val} = $vals[0];
	} elsif (@headers == 2) {
	    $input{$sheet}{$ref_val}{$vals[0]} = $vals[1];
	} elsif (@headers == 3) {
	    $input{$sheet}{$ref_val}{$vals[0]}{$headers[-2]} = $vals[-2];
	    $input{$sheet}{$ref_val}{$vals[0]}{$headers[-1]} = $vals[-1];
	} else {
	    die "Can't deal with more than 4 values!";
	}
    }
}


# combine collected values for output

# set up header line
my @header = ('#', 'Authors', 'Presenter', 'Title', 'Decision', 'Weighted score', 'Average quality score', 'Quality');

# score entries and mean average:
my @score_fields = ();
foreach my $i (sort { $a <=> $b } keys %{$input{'Review form fields'}}) {
    foreach my $field (keys %{$input{'Review form fields'}{$i}}) {
	my $add = $input{'Review form fields'}{$i}{$field}{'has scores?'};
	if ($add) {
	    # skip Quality scores - already listed at the start
	    next if ($field eq 'Quality score');
	    push @header, "$field scores", "$field mean";
	    push @score_fields, $i;
	}
    }
}

# comment entries
my @comment_fields = ();
foreach my $i (sort { $a <=> $b } keys %{$input{'Review form fields'}}) {
    foreach my $field (keys %{$input{'Review form fields'}{$i}}) {
	my $add = $input{'Review form fields'}{$i}{$field}{'has text?'};
	if ($add) {
	    push @header, "$field remarks";
	    push @comment_fields, $i;
	}
    }
}

# additional columns
push @header, ('Demo?', 'Formats', 'Flipped', 'Timezone', 'Diversity', 'Comments (Nomi)', 'Time (est.)', 'Comments from others');

# generate output as XLSX file:
my $col = 0;
my $row = 0;

# print header
foreach (@header) {
    $worksheet->write( $row, $col, $_, $formatb );
    $col++;
}

# one row per submission
foreach my $id (sort { $a <=> $b } keys %{$input{'Submissions'}}) {
    #    my $title = $input{'Submissions'}{$id};
    my @keys = sort keys %{$input{'Submissions'}{$id}};
    my $title = $keys[0];
    my $authors = $input{'Submissions'}{$id}{$title}{'authors'};
    my $presenter = $input{'Submission field values'}{$id}{'Presenter'};
    my $decision = $input{'Submissions'}{$id}{$title}{'decision'};
    my @out = ($id, $authors, $presenter, $title, $decision, '');
    
    # total scores (same as quality) and mean average
    my %review_ids = ();
    if (defined $input{'Reviews'}{$id}) {
	my @scores = ();
	foreach my $review (sort { $a <=> $b } keys %{$input{'Reviews'}{$id}}) {
	    push @scores, $input{'Reviews'}{$id}{$review};
	    $review_ids{$review}++;
	}
	my $scores = join ',', @scores;
	my $sum = 0;
	foreach (@scores) {
	    $sum += $_;
	}
	my $mean = sprintf "%.1f", $sum / @scores;
	push @out, $mean, $scores;
    } else {
	push @out, '', '';
    }

    # scores and mean average per individual submission field
    my @reviews = sort { $a <=> $b } keys %review_ids;
    foreach my $i (@score_fields) {
	my @scores = ();
	my $sum = 0;
	foreach my $review (@reviews) {
	    if (defined $input{'Review scores'}{$review}
		and
		defined $input{'Review scores'}{$review}{$i}) {
		push @scores, $input{'Review scores'}{$review}{$i};
		$sum += $input{'Review scores'}{$review}{$i};
	    }
	}
	my $scores = join ',', @scores;
	my $mean = '';
	if (@scores) {
	    $mean = sprintf "%.1f", $sum / @scores;
	}
	push @out, $scores, $mean;
    }

    my $wrapped = @out;
    
    # reviewers' comments per field
    foreach my $i (@comment_fields) {
	my @comments = ();
	foreach my $review (@reviews) {
	    if (defined $input{'Review texts'}{$review}
		and
		defined $input{'Review texts'}{$review}{$i}) {
		if ($input{'Review texts'}{$review}{$i}) {
		    push @comments, "$review: " . $input{'Review texts'}{$review}{$i};
		}
	    }
	}
	my $comments = join "\n", @comments;
	push @out, $comments;
    }

    # requested presentation format
    my $format = '';
    my $demo = '';
    my $flipped = '';
    my $timezone = '';
    my $diversity = '';
    if (defined $input{'Submission field values'}{$id}) {
	if (defined $input{'Submission field values'}{$id}{'Type'}) {
	    $format = $input{'Submission field values'}{$id}{'Type'};
	    $flipped = $input{'Submission field values'}{$id}{'flipped'};
	    $timezone = $input{'Submission field values'}{$id}{'Timezone'};
	    $diversity = $input{'Submission field values'}{$id}{'diversity'};
	    $format =~ s/, /\n/g;
	    if ($format =~ /demo/) {
		$demo = 'y';
	    }
	}
    }
    push @out, $demo, $format, $flipped, $timezone, $diversity;

    
    # write each cell to file:
    $row++;
    my $col = 0;
    foreach (@out) {
	my $format = '';

	# center scores
	if ($col >= 2) {
	    $format = $formatc;
	}

	# comments could be multi-line output, allow text wrapping
	if ($col >= $wrapped) {
            $format = $formatw;
	}

	# translate HTML code into symbols
	if (defined $_) {
	    s/\&lt;/</g;
	    s/\&gt;/>/g;
	    s/\&amp;/&/g;
	    s/\&apos;/'/g;
	    s/\&quot;/"/g;
	}
	$worksheet->write( $row, $col, $_, $format);
	$col++;
    }
}

$workbook->close();

warn "Output written to $out_file\n";


