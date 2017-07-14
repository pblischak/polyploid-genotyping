#!/usr/bin/perl -w
use strict;

open TOT_FILE, ">", $ARGV[0] || die "$!";
open REF_FILE, ">", $ARGV[1] || die "$!";

my $allele_depth_col = $ARGV[2];
my $depth_filter     = $ARGV[3];
my @lines = <STDIN>;

foreach my $line (@lines){
  chomp($line);
  if ($line =~ /^\s*#/) { next; }

  my @pieces = split /\t/, $line;
  my $pieces_size = @pieces;

  for(my $a = 9; $a < $pieces_size; $a = $a + 1){
    my $record = $pieces[$a];
    if($record eq "./."){
      print TOT_FILE "-9\t";
      print REF_FILE "-9\t";
      next;
    }
    my @record_entries = split /:/, $record;
    my $AD = $record_entries[$allele_depth_col - 1];
    my @AD_vals = split /,/, $AD;

    if($AD_vals[0] eq "."){
      print TOT_FILE "-9\t";
      print REF_FILE "-9\t";
    } else {
      my $AD_sum = $AD_vals[0] + $AD_vals[1];
      if($AD_sum < $depth_filter){
        print TOT_FILE "-9\t";
        print REF_FILE "-9\t";
      } else {
        print TOT_FILE "$AD_sum\t";
        print REF_FILE "$AD_vals[1]\t";
      }
    }
  }

  print TOT_FILE "\n";
  print REF_FILE "\n";
}
