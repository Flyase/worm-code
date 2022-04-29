#==============================================================
# File Name: bowtie2_depth.pl
# Author: wangyifeng
# mail: wangyifeng
# Created Time: Wed 25 Mar 2020 12:16:39 PM CST
#=============================================================
#!/usr/bin/perl -w
use strict;
die "perl $0 <genome> <_1.fastq.gz> <_2.fastq.gz>" unless @ARGV==3;

my $g = $ARGV[0];
my $r1 = $ARGV[1];
my $r2 = $ARGV[2];
#my $base = basename($r1);

print "date\n";
print "bowtie2-build $g $g\n";
print "bowtie2 --very-sensitive-local -x $g -1 $r1 -2 $r2 -p 24 | samtools view -b -@ 24 -q 33 - -o $r1.bam\n";
print "samtools sort -@ 24 -m 5g -o $r1.sort.bam $r1.bam\n";
print "samtools faidx $g\n";
print "cut -f 1,2 $g.fai \> $g.fai.g\n";
print "bedtools makewindows -g $g.fai.g -w 10000 | awk \'{print \$1\"\\t\"\$2\"\\t\"\$3}\' | sort -k1,1 -k2,2n \> $g.10k\n";
#samtools depth -m 180 -Q 59
print "samtools depth $r1.sort.bam | awk \'{print \$1\"\\t\"\$2-1\"\\t\"\$2\"\\t\"\$3}\' | sort -k1,1 -k2,2n | bedtools map -a  $g.10k -b - -c 4 -o sum,count,median,mean \> $r1.depth.10k-win\n";
print "date\n";
