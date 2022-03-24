#==============================================================
## File Name: htseq_count.pl
## Author: wangyifeng
## mail: wangyifeng
## Created Time: Tue 18 Aug 2020 12:16:39 PM CST
##=============================================================
##!/usr/bin/perl -w
#use strict;
### ',\,$,",>都需要\转义符，{}不需要
#htseq-count命令的常用参数： -r参数指定文件的排序方式，pos:按照染色体位置排序，name:按照read名称进行排序。双端测序数据必选参数，默认值为name。对于单端测序数据，该选项可以忽略。对于双端测序数据，必须要对SAM文件进行排序，对read name或 染色体位置 进行排序皆可，但是推荐使用read name进行排序，将会大大节约时间及服务器资源。通过 -r 可以指定数据是以什么方式排序的： 可以是 name 或 pos ， 默认是name.
####htseq-count -f bam -i transcript_id/Parent
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

