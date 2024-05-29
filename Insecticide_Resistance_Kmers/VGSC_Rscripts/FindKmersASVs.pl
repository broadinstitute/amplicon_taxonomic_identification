#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw($opt_k $opt_f $opt_o $opt_r);

# Usage
my $usage = "
FindKmersASVs.pl - Searches for a sets of kmers in a set of fasta sequences
Copyright (C) 2023 by Jacob A Tennessen 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Usage: perl FindKmersASVs.pl options
 required:
  -k  (path to) a kmer file [two columns, first column is kmer sequence, second column is category (e.g. 'Resistant')]
  -f  (path to) a sequence file in FASTA format
  -o  (path to) output file
optional:
  -r  ignore reverse complement
";

#############

# command line processing.
getopts('k:f:o:r');
die $usage unless ($opt_k);
die $usage unless ($opt_f);
die $usage unless ($opt_o);

my ($kmerfile, $seqfile, $outfile, $ignorerev);

$kmerfile = $opt_k if $opt_k;
$seqfile = $opt_f if $opt_f;
$outfile = $opt_o if $opt_o;

if (defined $opt_r) {
  $ignorerev = 1;
}

my %kmers;

my %counts;

open(FIRST, $kmerfile) || die "can't open $kmerfile\n";

while (<FIRST>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    if ($line =~ /\t/) {
      my @data = split "\t", $line;
      my $kmer = uc($data[0]);
      if ((defined $data[1])&&(($kmer =~ /A/)||($kmer =~ /C/)||($kmer =~ /G/)||($kmer =~ /T/))) {
        $kmers{$kmer} = $data[1];
        $counts{$data[1]} = 0;
      }
    }
}

close (FIRST);

my @nucs;

my @matches;

push @matches, "ASV\tCategory";

my $name;

open(SECOND, $seqfile) || die "can't open $seqfile\n";

while (<SECOND>) {
  my $line = $_;
  $line =~ s/\r|\n//g;
  if ($line =~ />/) {
    if ((defined $name)&&(defined $nucs[0])) {
      my $seq = uc(join "", @nucs);
      my $revcomp = reverse $seq;
      $revcomp =~ tr/ATGCatgc/TACGtacg/;
      my @foundkmers;
      foreach my $k (keys %kmers) {
        if (defined $ignorerev) {
          if ($seq =~ /$k/) {
            $counts{$kmers{$k}} += 1;
            push @foundkmers, $k;
          }
        } else {
          if (($seq =~ /$k/)||($revcomp =~ /$k/)) {
            $counts{$kmers{$k}} += 1;
            push @foundkmers, $k;
          }
        }
      }
      my $match = "No kmer found";
      if (defined $foundkmers[0]) {
        if (defined $foundkmers[1]) {
          print "WARNING: multiple kmers found in the same sequence. Kmers do not uniquely identify target region.\n";
          $match = "Multiple";
        } else {
          $match = $kmers{$foundkmers[0]};
        }
      }
      push @matches, "$name\t$match";
    }
    @nucs = ();
    my @namedata = split ">", $line;
    $name = $namedata[1];
  } else {
    push @nucs, $line;
  }
}

close (SECOND);

if ((defined $name)&&(defined $nucs[0])) {
  my $seq = uc(join "", @nucs);
  my $revcomp = reverse $seq;
  $revcomp =~ tr/ATGCatgc/TACGtacg/;
  my @foundkmers;
  foreach my $k (keys %kmers) {
    if (defined $ignorerev) {
      if ($seq =~ /$k/) {
        $counts{$kmers{$k}} += 1;
        push @foundkmers, $k;
      }
    } else {
      if (($seq =~ /$k/)||($revcomp =~ /$k/)) {
        $counts{$kmers{$k}} += 1;
        push @foundkmers, $k;
      }
    }
  }
  my $match = "No kmer found";
  if (defined $foundkmers[0]) {
    if (defined $foundkmers[1]) {
      print "WARNING: multiple kmers found in the same sequence. Kmers do not uniquely identify target region.\n";
      $match = "Multiple";
    } else {
      $match = $kmers{$foundkmers[0]};
    }
  }
  push @matches, "$name\t$match";
}

print "Category\tCount\n";

foreach my $c (sort keys %counts) {
  print "$c\t$counts{$c}\n";
}

my $result = join "\n", @matches;

unless ( open(OUT, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}
print OUT $result;
close (OUT);



