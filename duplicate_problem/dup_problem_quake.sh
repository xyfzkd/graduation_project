path0=/home/xieyufeng/test;
########QC1########
mkdir $path0/01.trim
mkdir $path0/01.trim/QC1
for i in `ls $path0/00.rawdata`;
do {j}=${i%.*};
echo ${j} fastqc;
fastqc ${j};
done
mv *.{html,zip} $path0/01.trim/QC1


########cutadpater########
path0=/home/xieyufeng/test;
cd $path0/02.mapping;
#rm -rf *;
cd $path0/01.trim
mkdir short cutadapt
for i in `ls $path0/00.rawdata`;
do j=${i%??.*};
echo ${j} >> $path0/filelist.txt
done
uniq $path0/filelist.txt $path0/list.txt

for j in `cat $path0/list.txt`;do
echo ${j} for_cutadapt;
cutadapt -j 10 -u 3 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -a AAAAAAAA -o $path0/01.trim/cutadapt/${j}_1.cutadapt.fastq $path0/00.rawdata/${j}_1.fastq
echo ${j} cut_6N;
cutadapt -j 10 -u -6 --trim-n -o $path0/01.trim/cutadapt/${j}_1.cut_6N.fastq $path0/01.trim/cutadapt/${j}_1.cutadapt.fastq
echo ${j} re_cutadapt;
cutadapt -j 10 -u 6 -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -a TTTTTTTT -o $path0/01.trim/cutadapt/${j}_2.cutadapt.fastq $path0/00.rawdata/${j}_2.fastq
echo ${j} cut_3N;
cutadapt -j 10 -u -3 --trim-n -o $path0/01.trim/cutadapt/${j}_2.cut_3N.fastq $path0/01.trim/cutadapt/${j}_2.cutadapt.fastq
echo ${j} filter;
cutadapt -j 10 -a GGGGGGGGGGGGGGGGGGGG -A GGGGGGGGGGGGGGGGGGGG -m 20 --too-short-paired-output=$path0/test/01.trim/short/$j.tooShort.fastq -o $path0/01.trim/cutadapt/${j}_1.filter.fastq -p $path0/01.trim/cutadapt/${j}_2.filter.fastq $path0/01.trim/cutadapt/${j}_1.cut_6N.fastq $path0/01.trim/cutadapt/${j}_2.cut_3N.fastq
done

########rm rRNA########
mkdir $path0/02.mapping/no_rRNA
cd $path0/02.mapping/no_rRNA
mkdir fastq sam rsem_bam
for j in `cat $path0/list.txt`;do
echo $j remove rRNA;
bowtie2 -p 20 --norc --sensitive-local --no-unal --un-conc $path0/02.mapping/no_rRNA/fastq/$j.no_rRNA.fq -x $path0/rRNA_index/rRNA -1 $path0/01.trim/cutadapt/${j}_1.filter.fastq -2 $path0/01.trim/cutadapt/${j}_2.filter.fastq -S $path0/02.mapping/no_rRNA/sam/$j.rRNA.sam
#rsem-tbam2gbam $path0/rRNA_index/rRNA $path0/02.mapping/no_rRNA/sam/$j.rRNA.sam $path0/02.mapping/no_rRNA/rsem_bam/$j.rRNA.rsem.clean.bam
done

########rm miRNA piRNA Y_RNA srpRNA tRNA snRNA snoRNA lncRNA mRNA tucp########
mkdir /home/xieyufeng/test/02.mapping/sequential_mapping
arr=("miRNA" "piRNA" "Y_RNA" "srpRNA" "tRNA" "snRNA" "snoRNA" "lncRNA" "mRNA" "tucp")
for k in {0..9};do
mkdir $path0/02.mapping/sequential_mapping/no_${arr[$k]}
cd $path0/02.mapping/sequential_mapping/no_${arr[$k]}
mkdir fastq sam rsem_bam
for j in `cat $path0/list.txt`;do
echo ${arr[$k]} mapping
let "l=k-1"
if [ $l==-1 ];then
bowtie2 -p 20 --sensitive-local --norc --no-unal --un-conc $path0/02.mapping/sequential_mapping/no_${arr[$k]}/fastq/$j.${arr[$k]}.unAligned.fastq -x $path0/RNA_index/${arr[$k]} -1 $path0/02.mapping/no_rRNA/fastq/$j.no_rRNA.1.fq -2 $path0/02.mapping/no_rRNA/fastq/$j.no_rRNA.2.fq -S $path0/02.mapping/sequential_mapping/no_${arr[$k]}/sam/$j.${arr[$k]}.sam
rsem-tbam2gbam -p 20 $path0/RNA_index/${arr[$k]} $path0/02.mapping/sequential_mapping/no_${arr[$k]}/sam/$j.${arr[$k]}.sam $path0/02.mapping/sequential_mapping/no_${arr[$k]}/rsem_bam/$j.${arr[$k]}.rsem.clean.bam 
else
bowtie2 -p 20 --sensitive-local --norc --no-unal --un-conc $path0/02.mapping/sequential_mapping/no_${arr[$k]}/fastq/$j.${arr[$k]}.unAligned.fastq -x $path0/RNA_index/${arr[$k]} -1 $path0/02.mapping/sequential_mapping/no_${arr[$l]}/fastq/$j.${arr[$l]}.unAligned.1.fastq -2 $path0/02.mapping/sequential_mapping/no_${arr[$l]}/fastq/$j.${arr[$l]}.unAligned.2.fastq -S $path0/02.mapping/sequential_mapping/no_${arr[$k]}/sam/$j.${arr[$k]}.sam
rsem-tbam2gbam -p 20 $path0/RNA_index/${arr[$k]} $path0/02.mapping/sequential_mapping/no_${arr[$k]}/sam/$j.${arr[$k]}.sam $path0/02.mapping/sequential_mapping/no_${arr[$k]}/rsem_bam/$j.${arr[$k]}.rsem.clean.bam 
fi
done
done

arr=("miRNA" "piRNA" "Y_RNA" "srpRNA" "tRNA" "snRNA" "snoRNA" "lncRNA" "mRNA" "tucp")
for k in {0..9};do
for j in `cat $path0/list.txt`;do
echo ${arr[$k]} tbam2gbam
#rsem-tbam2gbam $path0/RNA_index/${arr[$k]} $path0/02.mapping/sequential_mapping/no_${arr[$k]}/sam/$j.${arr[$k]}.sam $path0/02.mapping/sequential_mapping/no_${arr[$k]}/rsem_bam/$j.${arr[$k]}.rsem.clean.bam 
samtools view -@ 3 -S -b $path0/02.mapping/sequential_mapping/no_${arr[$k]}/sam/$j.${arr[$k]}.sam > $path0/02.mapping/sequential_mapping/no_${arr[$k]}/sam/$j.${arr[$k]}.bam
done
done

########mapping########
mkdir $path0/02.mapping/STAR_mapping
for j in `cat $path0/list.txt`;do
echo $j STAR mapping;
STAR --genomeDir $path0/STAR_index/ --runThreadN 20 --outFilterMismatchNoverLmax 0.05 --readFilesIn $path0/02.mapping/no_rRNA/fastq/$j.no_rRNA.1.fq $path0/02.mapping/no_rRNA/fastq/$j.no_rRNA.2.fq --outFileNamePrefix $path0/02.mapping/STAR_mapping/ --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outReadsUnmapped Fastx
done

cd $path0/02.mapping/
mkdir $path0/02.mapping/bowtie2_mapping
mkdir $path0/02.mapping/bowtie2_mapping/sam $path0/02.mapping/bowtie2_mapping/unmapped $path0/02.mapping/bowtie2_mapping/bam
for j in `cat $path0/list.txt`;do
echo $j bowtie2 mapping;
bowtie2 -p 20 --sensitive-local --no-unal --un-conc $path0/02.mapping/bowtie2_mapping/unmapped/$j.unmapped.fastq -x $path0/Bowtie2_index/GRCh38_p10 -1 $path0/02.mapping/no_rRNA/fastq/$j.no_rRNA.1.fq -2 $path0/02.mapping/no_rRNA/fastq/$j.no_rRNA.2.fq  -S $path0/02.mapping/bowtie2_mapping/sam/$j.bowtie2.sam > $path0/02.mapping/bowtie2_mapping/{bowtie2_mapping.log}
samtools view -@ 3 -S -b $path0/02.mapping/bowtie2_mapping/sam/$j.bowtie2.sam > $path0/02.mapping/bowtie2_mapping/bam/$j.bowtie2.bam
done
###MarkDuplicate###
for j in `cat $path0/list.txt`;do
echo $j bam sorting;
samtools sort -@ 8 $path0/02.mapping/bowtie2_mapping/bam/$j.bowtie2.bam > $path0/02.mapping/bowtie2_mapping/bam/$j.bowtie2.sorted.bam
done

###MarkDuplicate###
mkdir $path0/02.mapping/bowtie2_mapping/dedup
cd $path0/02.mapping/bowtie2_mapping/dedup
mkdir mark dedup metrics
for j in `cat $path0/list.txt`;do
echo $j
java -jar $path0/picard.jar MarkDuplicates I=$path0/02.mapping/bowtie2_mapping/bam/$j.bowtie2.sorted.bam  O=$path0/02.mapping/bowtie2_mapping/dedup/mark/$j.bowtie2.marked.bam M=$path0/02.mapping/bowtie2_mapping/dedup/metrics/$j.marked.metrics.txt
done

for j in `cat $path0/list.txt`;do
echo $j
java -jar $path0/picard.jar MarkDuplicates I=$path0/02.mapping/bowtie2_mapping/bam/$j.bowtie2.sorted.bam  O=$path0/02.mapping/bowtie2_mapping/dedup/dedup/$j.bowtie2.dedup.bam M=$path0/02.mapping/bowtie2_mapping/dedup/metrics/$j.marked.dedup.txt REMOVE_SEQUENCING_DUPLICATES=true
done
###Make table###
mkdir $path0/02.mapping/bowtie2_mapping/dedup
mkdir $path0/02.mapping/bowtie2_mapping/table
for j in `cat $path0/list.txt`;do
echo $j
samtools view $path0/02.mapping/bowtie2_mapping/dedup/mark/$j.bowtie2.marked.bam | cut -f 1,2,3,4 > $path0/02.mapping/bowtie2_mapping/table/$j.csv
done

###Read type###
for i in "miRNA" "piRNA" "Y_RNA" "srpRNA" "tRNA" "snRNA" "snoRNA" "lncRNA" "mRNA" "tucp";do
perl $path0/bed/gtf2Bed.pl $path0/GTF/$i.gtf > $path0/bed/$i.bed
done

mkdir $path0/02.mapping/bowtie2_mapping/type
for i in "miRNA" "piRNA" "Y_RNA" "srpRNA" "tRNA" "snRNA" "snoRNA" "lncRNA" "mRNA" "tucp";do
for j in `cat $path0/list.txt`;do
echo $j $i anno
samtools view -bL $path0/bed/$i.bed $path0/02.mapping/bowtie2_mapping/dedup/mark/$j.bowtie2.marked.bam > $path0/02.mapping/bowtie2_mapping/type/$j.$i.bam
done
done

mkdir $path0/02.mapping/bowtie2_mapping/table_read
for j in `cat $path0/list.txt`;do
for i in "miRNA" "piRNA" "Y_RNA" "srpRNA" "tRNA" "snRNA" "snoRNA" "lncRNA" "mRNA" "tucp";do
echo $j $i cut
samtools view -@ 4 $path0/02.mapping/bowtie2_mapping/type/$j.$i.bam | cut -f 1,2,3,4 > $path0/02.mapping/bowtie2_mapping/table_read/$j.$i.csv
done
done

for i in "miRNA" "piRNA" "Y_RNA" "srpRNA" "tRNA" "snRNA" "snoRNA" "lncRNA" "mRNA" "tucp";do
perl $path0/bed/gtf2Bed.pl $path0/GTF/$i.gtf >> $path0/bed/all.bed
done

for j in `cat $path0/list.txt`;do
echo $j all anno
samtools view -bL $path0/bed/all.bed $path0/02.mapping/bowtie2_mapping/dedup/mark/$j.bowtie2.marked.bam > $path0/02.mapping/bowtie2_mapping/type/$j.bam
done 

for j in `cat $path0/list.txt`;do
echo $j sort
samtools sort -@ 8 $path0/02.mapping/bowtie2_mapping/type/$j.bam > $path0/02.mapping/bowtie2_mapping/type/$j.sorted.bam 
done

for j in `cat $path0/list.txt`;do
echo $j anno
samtools view -@ 4 $path0/02.mapping/bowtie2_mapping/type/$j.sorted.bam | cut -f 1,2,3,4 > $path0/02.mapping/bowtie2_mapping/table_read/$j.csv
done
