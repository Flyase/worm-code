samtools faidx Tmur.fa
/usr/local/jdk1.8.0_101/bin/java -jar /share/app/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R= Tmur.fa O= Tmur.dict
