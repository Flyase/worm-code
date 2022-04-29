#==============================================================
# File Name: snpcalling.pl
# Author: wangyifeng
# mail: wangyifeng
# Created Time: Tue 07 Apr 2020 12:47:08 PM CST
#=============================================================
#!/usr/bin/perl -w
use strict;
use File::Basename;

die "perl $0 <genome> <R1.gz> <R2.gz>" unless @ARGV==3;


my $g = $ARGV[0];
my $r1 = $ARGV[1];
my $r2 = $ARGV[2];
my $base = basename($r1);

print "### bwa-map,best practices mapping using bwa+picard\n";
print "date\n";
print "bwa index $g -p $g\n";
print "bwa mem -t 14 -M $g $r1 $r2 | samtools view -bS - > $r1.bam\n";
print "date\n";
print "### depth\n";
print "samtools sort -@ 24 -m 5g -o $r1.sort.bam $r1.bam\n";
print "samtools faidx $g\n";
print "cut -f 1,2 $g.fai \> $g.fai.g\n";
print "bedtools makewindows -g $g.fai.g -w 50000 -s 40000 | awk \'{print \$1\"\\t\"\$2\"\\t\"\$3}\' | sort -k1,1 -k2,2n \> $g.50k\n";
print "samtools depth $r1.sort.bam | awk \'{print \$1\"\\t\"\$2-1\"\\t\"\$2\"\\t\"\$3}\' | sort -k1,1 -k2,2n | bedtools map -a  $g.50k -b - -c 4 -o sum,count,median,mean \> $r1.depth.50k-win\n";
print "date\n";
print "### SNP-1\n";
print "samtools view -bS $r1.bam -t $g.fai -o $r1.fai.bam\n";
print "samtools sort --threads 14 -m 2000000000 $r1.fai.bam -o $r1.fai.sort.bam\n";
print "java -Xmx50G -jar picard.jar AddOrReplaceReadGroups I=$r1.fai.sort.bam O=$r1.fai.sort.group.bam RGLB=$base RGPL=illumina RGSM=$base RGPU=$base\n";
print "java -Xmx50G -jar picard.jar SortSam I=$r1.fai.sort.group.bam O=$r1.fai.sort.group.sort.bam SORT_ORDER=coordinate\n";
print "/java -Xmx50G -jar picard.jar MarkDuplicates I=$r1.fai.sort.group.sort.bam O=$r1.fai.sort.group.sort.dup.bam METRICS_FILE=$r1.metrics.txt MAX_FILE_HANDLES=800\n";
print "samtools index $r1.fai.sort.group.sort.dup.bam\n";
print "date\n";
print "### SNP-2,indel realignment,BaseRecalibrator\n";
print "java -Xmx50G -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R $g -I $r1.fai.sort.group.sort.dup.bam -o $r1.fai.sort.group.sort.dup.bam.intervals 1\>intervals2.log 2\>\&1\n";
print "java -Xmx50G -jar GenomeAnalysisTK.jar -T IndelRealigner -R $g -I $r1.fai.sort.group.sort.dup.bam -targetIntervals $r1.fai.sort.group.sort.dup.bam.intervals -o $r1.fai.sort.group.sort.dup.realign.bam 1\>realign2.log 2\>\&1\n";
print "java -Xmx50G -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R $g --run_without_dbsnp_potentially_ruining_quality -I $r1.fai.sort.group.sort.dup.realign.bam  -o $r1.fai.sort.group.sort.dup.realign.recal.grp 1\>recal2.log 2\>\&1\n";
print "java -Xmx50G -jar GenomeAnalysisTK.jar -T PrintReads  -R $g -I $r1.fai.sort.group.sort.dup.realign.bam   -BQSR $r1.fai.sort.group.sort.dup.realign.recal.grp -o $r1.fai.sort.group.sort.dup.realign.recal.bam --num_bam_file_handles 100 1\>printRead2.log 2\>\&1\n";
print "date\n";
print "#HaplotypeCaller call snp\n";
print "java -Xmx50G -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R $g -ERC GVCF -I $r1.fai.sort.group.sort.dup.realign.recal.bam -o $r1.g.vcf.gz\n";
print "###genotype\n";
print "/public1/home/yifeng/yifeng/java/jre1.8.0_121/bin/java -Xmx50G -jar GenomeAnalysisTK.jar -R $g -T GenotypeGVCFs --variant $r1.g.vcf.gz -o $r1.g.raw.vcf.gz\n";
print "###filter snp\n";
print "perl /public1/home/yifeng/yifeng/code/pipeline/SNP/QualFilter.pl $r1.g.raw.vcf.gz $r1.g.QualFilter.vcf\n";
print "java -Xmx50G -jar GenomeAnalysisTK.jar -T SelectVariants -R $g -V $r1.g.QualFilter.vcf -selectType SNP -o $r1.snp.vcf\n";
print "java -Xmx50G -jar GenomeAnalysisTK.jar -T VariantFiltration -R $g -V $r1.snp.vcf --logging_level ERROR --filterExpression \" QD \< 2.0 || FS \> 60.0 || MQRankSum \< -12.5 || ReadPosRankSum \< -8.0 || SOR \> 3.0 || MQ \< 40.0 \" --filterName \"my_snp_filter\" -o $r1.snp.filter.vcf\n";
print "date\n";
