cd ./Data/Ref/hg19
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr2.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr3.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr4.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr5.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr6.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr7.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr8.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr9.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr10.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr11.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr12.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr13.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr14.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr15.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr16.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr17.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr18.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr20.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr21.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chrX.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chrY.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chrM.fa.gz
gunzip *gz
cat chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chr20.fa chr21.fa chr22.fa chrX.fa chrY.fa chrM.fa >hg19.fa
bowtie2-build hg19.fa hg19

cd ../mm9
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr1.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr2.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr3.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr4.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr5.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr6.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr7.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr8.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr9.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr10.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr11.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr12.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr13.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr14.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr15.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr16.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr17.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr18.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr19.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chrX.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chrY.fa.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chrM.fa.gz
gunzip *gz
cat chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chrX.fa chrY.fa chrM.fa >mm9.fa
bowtie2-build mm9.fa mm9
