#!/usr/bin/perl -w
die "perl $0 <reads.dir> <index> <out.dir>\nnote: reads must paired and merged\n" unless @ARGV==4; 
use Cwd qw(abs_path);
use File::Spec;
###read files
opendir DIR, $ARGV[0];
my @file = readdir DIR;
my $path = abs_path($ARGV[0]);
my $out_path = abs_path($ARGV[3]);
#my $index = abs_path($ARGV[1]); #linked file will point to its original path
my $index = File::Spec->rel2abs($ARGV[1]); #this one point to its absolute path
my $gtf = File::Spec->rel2abs($ARGV[2]); #this one point to its absolute path
###restore files
my %hash;
foreach my $file (@file) {
        next unless $file =~ /(.*)\.(fq|fastq)\.gz$/;
        my $read_name = $1;
        $hash{$read_name}=1;
};
###write work.sh
open OUT, ">$out_path/run_htcout.sh";
foreach my $read_name (sort keys %hash) {
        print OUT "hisat2 -x $index -U $path/$read_name.fastq.gz -p 12 | samtools view -q 33 -b -@ 12 - -o $out_path/$read_name.bam\n"; #samtools view -u，samtools sort -
        print OUT "samtools sort -n -@ 14 -m 5g -o $out_path/$read_name.sort.bam $out_path/$read_name.bam\n";
        print OUT "htseq-count -f bam -r name -i transcript_id -t exon $out_path/$read_name.sort.bam $gtf -s no > $out_path/$read_name.rno.count\n"; 
};
close OUT;
