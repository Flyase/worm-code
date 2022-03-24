#!/usr/bin/perl
use strict;

open(IN,"gzip -dc $ARGV[0] |") || die "Can't open IN!\n";
open(OUT,">$ARGV[1]") || die "Can't open OUT!\n";

while(<IN>){
	chomp;
	if(/^#/){
		print OUT "$_\n";
	}
	else{
		my @tmp=split /\s+/,$_;
		next if($tmp[5]<30);
		print OUT "$_\n";
	}
}
