# dup_problem
## for pico data
path0 = `/home/xieyufeng/test1`
## for quake data
path0 = `/home/xieyufeng/test`

## script
`dup_problem_quake.sh`

```
#parameter attention
bowtie2 --sensitive-local
#gtf2bed
perl $path0/bed/gtf2Bed.pl $path0/GTF/$i.gtf >> $path0/bed/all.bed
#overlapping to divide the bam file into different type bam files
samtools view -bL $path0/bed/all.bed $path0/02.mapping/bowtie2_mapping/dedup/mark/$j.bowtie2.marked.bam > $path0/02.mapping/bowtie2_mapping/type/$j.bam

```
## stat file
```
#csv files transformed from bam files, deleted some columns and only left chr、loci、1024、ID information
$path0/02.mapping/bowtie2_mapping/table_read/
```

## visulization
* for overall dispaly
* for different RNA type dispaly

`dup_read.ipynb`
other scripts written by xupeng, maybe in `home/chenxupeng/projects/exseek/jupyter/`
