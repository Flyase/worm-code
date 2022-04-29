###mkdir 01.map
cd 01.map
tophat -p 20 ../00.data/Bmal.fe.fa ../00.data/SRR4308249_1.fastq.gz,../00.data/SRR4308257_1.fastq.gz,../00.data/SRR4308252_1.fastq.gz ../00.data/SRR4308249_2.fastq.gz,../00.data/SRR4308257_2.fastq.gz,../00.data/SRR4308252_2.fastq.gz \
cd tophat_out \
samtools sort -n unmapped.bam -o unmapped.sort.bam \
bamToFastq -i unmapped.sort.bam -fq unmapped_1.fq -fq2 unmapped_2.fq 
###mkdir 02.trinity
cd 02.trinity

trinityrnaseq-Trinity-v2.4.0/Trinity \
	--seqType fq \
	-max_memory 20G \
	--left unmapped_1.fq \
	--right unmapped_2.fq \
	--CPU 20 \
	--output trinity_out_dir_all \

### ORF annotation
TransDecoder.LongOrfs -t Trinity.fasta -m 200
TransDecoder.Predict -t Trinity.fasta --single_best_only

###mkdir 03.align
cd 03.align/
ln -s ../02.trinity/trinity_out_dir_all/Trinity.fasta
ln -s ../00.data/Bmal.fe.fa
nohup nucmer --mum --prefix tri Bmal.fe.fa Trinity.fasta &
show-coords tri.delta > tri.delta.coords
bioawk -c fastx '{ print $name, length($seq) }' Trinity.fasta | sort -k1,1 > Trinity.fasta.length
##aligns
awk '{print $13,$4,$5,$10}' tri.delta.coords | sed 1,5d | awk '$2>$3{print $1,$3,$2,$4}' > 1
awk '{print $13,$4,$5,$10}' tri.delta.coords | sed 1,5d | awk '$2<$3{print $1,$2,$3,$4}' > 2
cat 1 2 | sort -k1,1 -k2,2n | sed 's/ /\t/g' > tri.delta.coords.id.bed
########filtration
##2.ID>=50 & qry aligned length >=50% awk '$4>50{print $1,$2,$3}' tri.delta.coords.id.bed | sort -k1,1 -k2,2n | sed 's/ /\t/g' | bedtools merge -i - | awk '{print $1,$3-$2}' | awk '{sum[$1]+=$2}END{for(c in sum){print c,sum[c]}}' | sort -k1,1 > tri.delta.coords.bed.length  
join tri.delta.coords.bed.length Trinity.fasta.length | awk '{print $0,$2/$3}' | awk '$4>0.5' | awk '{print ">"$1}' | sort | uniq > tri.delta.coords.idmore50.rm.contig
awk '/^>/{print "";}{ printf"%s",/^>/ ? $0" ":$0}' Trinity.fasta  | awk 'NR==FNR {a[$1]} NR>FNR&&!($1 in a)' tri.delta.coords.idmore50.rm.contig - | awk '{print $1"\n"$NF}' | grep -v '^$' > Trinity.keep.fasta

###mkdir 04.filt
cd 04.filt
ln -s ../03.align/Trinity.keep.fasta
bowtie2-build Trinity.keep.fasta Trinity.keep.fasta    #####long time
conda activate py2
nohup tophat -p 10 -o SRR4308246 Trinity.keep.fasta SRR4308246_1.fastq.gz SRR4308246_2.fastq.gz &   ##female rna-seq reads
nohup tophat -p 10 -o SRR4308254 Trinity.keep.fasta SRR4308254_1.fastq.gz SRR4308254_2.fastq.gz &
nohup tophat -p 10 -o SRR4308259 Trinity.keep.fasta SRR4308259_1.fastq.gz SRR4308259_2.fastq.gz &

##mkdir 1.fe_rm.depth
ln -s ../tophat_out/accepted_hits.bam
samtools sort -@2 -m 2g -o accepted_hits.sort.bam accepted_hits.bam
samtools merge fe.sorted.bam ../04.filt/SRR4308259/accepted_hits.sort.bam ../accepted_hits.sort.bam ../04.filt/SRR4308246/accepted_hits.sort.bam
samtools depth -Q 30 fe.sorted.bam | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' | sort -k1,1 -k2,2n >  fe_depth_per_base.bed
##keep per_base_count>5 loci
awk '$4>5{print $1}' fe_depth_per_base.bed | sort | uniq -c | awk '{print $2,$1}' | sort -k1,1 > fe.cov.length 
bioawk -c fastx '{ print $name, length($seq) }' Trinity.keep.fasta | sort -k1,1 > Trinity.keep.fasta.length
join fe.cov.length Trinity.keep.fasta.length | awk '{print $0,$2/$3}' | awk '$4<0.5' | awk '{print ">"$1}' > Trinity.afterfecov.keep.contig
awk '/^>/{print "";}{ printf"%s",/^>/ ? $0" ":$0}' Trinity.keep.fasta | awk 'NR==FNR{a[$1]=$0;next}NR>FNR{if($1 in a)print $0}' Trinity.afterfecov.keep.contig - | awk '{print $1"\n"$2}' > Trinity.afterfecov.keep.fasta
##mkdir 2.m_depth
ln -s ../1.fe_rm.depth/Trinity.afterfecov.keep.fasta
hisat2-build Trinity.afterfecov.keep.fasta Trinity.afterfecov.keep.fasta
hisat2 -x Trinity.afterfecov.keep.fasta -1 ERR279690_1.fastq.gz -2 ERR279690_2.fastq.gz -p 12 | samtools view -b -@ 12 - -o m.bam
samtools sort -@2 -m 2g -o m.sort.bam m.bam
samtools depth m.sort.bam | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' | sort -k1,1 -k2,2n >  m_depth_per_base.bed
##keep per_base_count>5 loci
awk '$4>5{print $1}' m_depth_per_base.bed | sort | uniq -c | awk '{print $2,$1}' | sort -k1,1 > m.cov.length 
ln -s ../../03.align/Trinity.keep.fasta
bioawk -c fastx '{ print $name, length($seq) }' Trinity.afterfecov.keep.fasta | sort -k1,1 > Trinity.afterfecov.keep.fasta.length
join m.cov.length Trinity.afterfecov.keep.fasta.length | awk '{print $0,$2/$3}' | awk '$4>0.5' | awk '{print ">"$1}' > Trinity.afterfemcov.keep.contig
awk '/^>/{print "";}{ printf"%s",/^>/ ? $0" ":$0}' Trinity.afterfecov.keep.fasta | awk 'NR==FNR{a[$1]=$0;next}NR>FNR{if($1 in a)print $0}' Trinity.afterfemcov.keep.contig - | awk '{print $1"\n"$2}' > Trinity.final.fasta

####mkdir 05.depth
ln -s ../00.data/Bmal.fa
awk '/^>/{print "";}{ printf"%s",/^>/ ? $0" ":$0}' Bmal.fa | grep -v "Wolbachia" | grep -v "Chr" | awk '{print $1"\n"$2}' | grep -v '^$' > Bmal.rest.fa
ln -s ../04.filt/1.fe_rm.depth/Trinity.afterfecov.keep.fasta Trinity.final.fasta
nohup nucmer --mum --prefix tri Bmal.rest.fa Trinity.final.fasta &
show-coords tri.delta > tri.delta.coords
bioawk -c fastx '{ print $name, length($seq) }' Trinity.final.fasta | sort -k1,1 > Trinity.fasta.length
bioawk -c fastx '{ print $name, length($seq) }' Bmal.rest.fa | sort -k1,1 > Bmal.rest.fa.length

##aligns
awk '{print $13,$4,$5,$10}' tri.delta.coords | sed 1,5d | awk '$2>$3{print $1,$3,$2,$4}' > 1
awk '{print $13,$4,$5,$10}' tri.delta.coords | sed 1,5d | awk '$2<$3{print $1,$2,$3,$4}' > 2
cat 1 2 | sort -k1,1 -k2,2n | sed 's/ /\t/g' > tri.delta.coords.id.bed
########filtration
##2.ID>=50 & qry aligned length >=50%, mapped contig, rm the correspoding scaf 
awk '$4>50{print $1,$2,$3}' tri.delta.coords.id.bed | sort -k1,1 -k2,2n | sed 's/ /\t/g' | bedtools merge -i - | awk '{print $1,$3-$2}' | awk '{sum[$1]+=$2}END{for(c in sum){print c,sum[c]}}' | sort -k1,1 > tri.delta.coords.bed.length  
join tri.delta.coords.bed.length Trinity.fasta.length | awk '{print $0,$2/$3}' | awk '$4>0.5' | awk '{print $1}' | sort | uniq > tri.delta.coords.idmore50.map.contig
awk 'NR==FNR{a[$1]=$0;next}NR>FNR{if($13 in a)print $0}' tri.delta.coords.idmore50.map.contig tri.delta.coords  | awk '{print ">"$(NF-1)}' | sort | uniq > rm.contig
awk '/^>/{print "";}{ printf"%s",/^>/ ? $0" ":$0}' Bmal.fa | awk 'NR==FNR {a[$1]} NR>FNR&&!($1 in a)' rm.contig - | sort -k1,1 | awk '{print $1"\n"$2}' > Bmal.fa1
cat Bmal.fa1 Trinity.final.fasta | grep -v '^$' > Bmal.combine.fa
