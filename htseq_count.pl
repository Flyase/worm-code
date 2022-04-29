#==============================================================
## File Name: htseq_count.pl
## Author: wangyifeng
## mail: wangyifeng
## Created Time: Tue 18 Aug 2020 12:16:39 PM CST
##=============================================================
##!/usr/bin/perl -w
#use strict;

my $g = $ARGV[0];
my $r1 = $ARGV[1];
my $r2 = $ARGV[2];
my $gtf = $ARGV[3];
##my $base = basename($r1);
print "date\n";
print "#/public/home/lijing/software/hisat2-2.0.4/hisat2-build $g -p $g\n";
print "/public/home/lijing/software/hisat2-2.0.4/hisat2 -x $g -1 $r1 -2 $r2 -p 12 | samtools view -q 33 -b -@ 12 - -o $r1.bam\n";
print "samtools sort -n -@ 14 -m 5g -o $r1.sort.bam $r1.bam\n";
print "/public1/home/yifeng/yifeng/tools/miniconda2/bin/htseq-count -f bam -r name -i transcript_id -t CDS $r1.sort.bam $gtf -s no > $r1.count\n";
print "date\n";

