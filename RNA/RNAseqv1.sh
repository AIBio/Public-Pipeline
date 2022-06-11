#!/bin/bash


### ---------------
### Whole procedure
### ---------------
# 1. Build the working directory (required);
# 2. Fastqc + MultiQC: check the quality of sequencing (required);
# 3. Trimgalore + Fastqc + MultiQC: trim reads and QC (required);
# 4. STAR: align the reads to genome by STAR (required);
# 5. FeatureCounts: quantify gene and repeat locus expression levels (required);


### ---------------
### Running command
### ---------------
# nohup bash RNAseqv1.sh -c RNAseqv1.conf &


### --------------
### Parse augments
### --------------
Help()
{
cat >&2 <<EOF
FUNCTION:
        Custome pipeline to perform Bulk RNA-seq data analysis;

USAGE:
        $0 -c <configuration file> -h
        -c <*.conf file>      : configuration file that record all files needed by this pipeline;
        -h <help document>    : help document;

CONF FILE FORMAT:
        1. workding directory;
        2. full path of raw data;
        3. txt file that records sample names;
        4. data type: support PE and SE data;
        5. genomic sequence *.fa file;
        6. gene GTF file;
        7. adapter type: illumina/nextera;
EOF
   exit 1
}
while getopts "c:h" opt
do
   case "$opt" in
      c ) cf="$OPTARG" ;;
      h ) Help ;;
   esac
done
if [ -z "${cf}" ]; then
   echo "Some or all of the parameters are empty"
   Help
fi
echo "> ---------------------------- <"
# 1. workding directory
wd=`grep -v "^#" $cf | awk '{if(NR==1) print $0}'`;            echo -e "[1] working directory: ${wd}"
# 2. full path of raw data
rawdata=`grep -v "^#" $cf | awk '{if(NR==2) print $0}'`;       echo -e "[2] directory of raw fastq files: ${rawdata}"
# 3. txt file that records sample names
samples=`grep -v "^#" $cf | awk '{if(NR==3) print $0}'`;       echo -e "[3] sample file: ${samples}"
# 4. data type: support PE and SE data
dt=`grep -v "^#" $cf | awk '{if(NR==4) print $0}'`;            echo -e "[4] data type: ${dt}"
# 5. genomic sequence *.fa file
genome_fa=`grep -v "^#" $cf | awk '{if(NR==5) print $0}'`;     echo -e "[5] Genomic sequence file: ${genome_fa}"
# 6. gene GTF file
gene_gtf=`grep -v "^#" $cf | awk '{if(NR==6) print $0}'`;      echo -e "[6] gene GTF annotation file: ${gene_gtf}"
# 7. adapter type: illumina/nextera
adapter=`grep -v "^#" $cf | awk '{if(NR==7) print $0}'`;       echo -e "[7] adapter type: ${adapter}"
echo "> ---------------------------- <"


### ---------
### Functions
### ---------
# Create soft links for raw data
DataLink()
{
   indir=$1; sample=$2; outdir=$3
   [ ! -d ${outdir} ] && mkdir -p ${outdir}
   ls ${indir} | grep -f ${sample} - | while read file
   do
       format=`echo ${file} | awk -F'[.]' '{print $(NF-1)"."$NF}'`
       if [ "${format}" == "fastq.gz" ]; then
          prefix=`echo ${file} | sed 's,.fastq.gz,,g'`
       elif [ "${format}" == "fq.gz" ]; then
          prefix=`echo ${file} | sed 's,.fq.gz,,g'`
       fi
       ln -s ${indir}/${file} ${outdir}/${prefix}.fq.gz
   done
}
# Fastqc
QcReads()
{
   indir=$1; outdir=$2; logdir=$3
   [ ! -d ${outdir} ] && mkdir -p ${outdir}; [ ! -d ${logdir} ] && mkdir -p ${logdir}
   find ${indir} -name "*.fq.gz" | xargs -P 10 -I{} fastqc {} -o ${outdir} > ${logdir}/run.log 2>&1
}
# Multiqc
ParseQc()
{
   indir=$1; outdir=$2; logdir=$3
   [ ! -d ${outdir} ] && mkdir -p ${outdir}; [ ! -d ${logdir} ] && mkdir -p ${logdir}
   multiqc ${indir} -o ${outdir} > ${logdir}/run.log 2>&1
}
# Trimgalore
TrimReads()
{
   indir=$1; file=$2; outdir=$3; logdir=$4; qual=$5; len=$6; ada=$7
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ ${dt} == "PE" ]; then
     for i in `seq 1 2 $(cat ${file} | wc -l)`
     do
         read1=$(awk -v row=${i} '(NR == row){print $0}' ${file})
         r1pre=`echo ${read1} | xargs basename -s ".fq.gz"`
         read2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file})
         r2pre=`echo ${read1} | xargs basename -s ".fq.gz"`
         trim_galore --paired ${indir}/${read1} ${indir}/${read2} -o ${outdir} \
                     --quality ${qual} --max_n 4 --length ${len} \
                     --${ada} --cores 16 > ${logdir}/${r1pre}_${r2pre}.log 2>&1
     done
   elif [ ${dt} == "SE" ]; then
     while read file
     do
         prefix=`echo ${file} | xargs basename -s ".fq.gz"`
         trim_galore ${indir}/${file} -o ${outdir} --quality ${qual} --max_n 4 \
                     --length ${len} --${ada} --cores 16 > ${logdir}/${prefix}.log 2>&1
     done < ${file}
   else
     echo "Unrecognized data type!"; exit 1
   fi
}
# Star
StarBxWithGtf(){
   seq=$1; gtf=$2; outdir=$3; logdir=$4
   [ ! -d "${outdir}" ] && mkdir -p ${outdir} && [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   STAR --runThreadN 20 --runMode genomeGenerate --genomeDir ${outdir} \
        --genomeFastaFiles ${seq} --sjdbGTFfile ${gtf} > ${logdir}/star_build_index_with_gtf.log 2>&1
}
StarGene()
{
   indir=$1; file=$2; outdir=$3; logdir=$4; indexdir=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ ${dt} == "PE" ]; then
     for i in `seq 1 2 $(cat ${file} | wc -l)`
     do
         r1=$(awk -v row=${i} '(NR == row){print $0}' ${file})
         r1name=`echo ${r1} | xargs basename -s ".fq.gz"`
         r2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file})
         r2name=`echo ${r2} | xargs basename -s ".fq.gz"`
         STAR --runMode alignReads --runThreadN 16 --readFilesCommand zcat --outBAMsortingThreadN 10 \
              --outSAMattributes NH HI NM MD XS nM AS MD jM jI MC ch \
              --outFileNamePrefix ${outdir}/${r1name}_${r2name} --outSAMtype BAM SortedByCoordinate \
              --outSAMmultNmax -1 --outFilterMultimapNmax 10 --genomeDir ${indexdir} \
              --readFilesIn ${indir}/${r1},${indir}/${r2} > ${logdir}/${r1name}_${r2name}.log 2>&1
     done
   elif [ ${dt} == "SE" ]; then
     while read fq
     do
         prefix=`echo ${fq} | xargs basename -s ".fq.gz"`
         STAR --runMode alignReads --runThreadN 16 --readFilesCommand zcat --outBAMsortingThreadN 10 \
              --outSAMattributes NH HI NM MD XS nM AS MD jM jI MC ch \
              --outFileNamePrefix ${outdir}/${prefix} --outSAMtype BAM SortedByCoordinate \
              --outSAMmultNmax -1 --outFilterMultimapNmax 10 --genomeDir ${indexdir} \
              --readFilesIn ${indir}/${fq} > ${logdir}/${prefix}.log 2>&1
     done < ${file}
   else
     echo "Unrecognized data type!"; exit 1
   fi
}
# Subread
CountGene()
{
   indir=$1; gtf=$2; outdir=$3; logdir=$4
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   bamlist=$(ls ${indir}/*.bam | tr '\n' ' '); echo ${bamlist}
   featureCounts -a ${gtf} -o ${outdir}/all_samples_gene_count_id_matrix.txt \
                 -g 'gene_id' -T 12 ${bamlist} > ${logdir}/all_samples_gene_count_gene_id_matrix.log 2>&1
   featureCounts -a ${gtf} -o ${outdir}/all_samples_gene_count_name_matrix.txt \
                 -g 'gene_name' -T 12 ${bamlist} > ${logdir}/all_samples_gene_count_gene_name_matrix.log 2>&1
}


### ------------
### Run pipeline
### ------------
echo "1st step: Organize the workding directory (required)"
echo "Begin at $(date)"
cd ${wd}
mkdir logs results rawdata metadata scripts
DataLink ${rawdata} ${samples} ${wd}/rawdata
echo "Finish at $(date)"



echo "2nd step: Quality control of raw data (required)"
echo "Begin at $(date)"
ls ${wd}/rawdata/ | grep "fq.gz$" > ${wd}/metadata/raw_fq_file.txt
QcReads ${wd}/rawdata ${wd}/results/fastqc/raw ${wd}/logs/fastqc/raw
ParseQc ${wd}/results/fastqc/raw ${wd}/results/multiqc/raw ${wd}/logs/multiqc/raw
echo "Finish at $(date)"



echo "3rd step: Filter and trim the reads (required)"
echo "Begin at $(date)"
quality=20; length=30
TrimReads ${wd}/rawdata ${wd}/metadata/raw_fq_file.txt ${wd}/results/trimgalore ${wd}/logs/trimgalore ${quality} ${length} ${adapter}
QcReads ${wd}/results/trimgalore ${wd}/results/fastqc/trimgalore ${wd}/logs/fastqc/trimgalore
ParseQc ${wd}/results/fastqc/trimgalore ${wd}/results/multiqc/trimgalore ${wd}/logs/multiqc/trimgalore
if [ "${dt}" == "PE" ]; then
   ls ${wd}/results/trimgalore/ | grep "val" | grep "fq.gz$" > ${wd}/metadata/trim_fq_${dt}.txt
else
   ls ${wd}/results/trimgalore/ | grep "trimmed" | grep "fq.gz$" > ${wd}/metadata/trim_fq_${dt}.txt
fi
echo "Finish at $(date)"



echo "4th step: Mapping (required)"
echo "Begin at $(date)"
index_dir=${wd}/results/star/index
StarBxWithGtf ${genome_fa} ${gene_gtf} ${index_dir} ${wd}/logs/star/index
StarGene ${wd}/results/trimgalore ${wd}/metadata/trim_fq_${dt}.txt ${wd}/results/star/gene ${wd}/logs/star/gene ${index_dir}
ParseQc ${wd}/results/star/gene ${wd}/results/multiqc/star_gene ${wd}/logs/multiqc/star_gene gene_mapping
echo "Finish at $(date)"



echo "5th step: Quantification (required)"
echo "Begin at $(date)"
CountGene ${wd}/results/star/gene ${gene_gtf} ${wd}/results/featurecounts/gene ${wd}/logs/featurecounts/gene
echo "Finish at $(date)"
