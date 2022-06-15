#!/bin/sh

### ---------------
### Whole procedure
### ---------------

# 1. Build the workding directories
# 2. Fastqc + MultiQC: check the quality of sequencing
# 3. Trimgalore + Fastqc + MultiQC: trim reads and QC
# 4. Bowtie2: create genome index and align the reads to genome by Bowtie2
# 5. Sambamba: re-format sam files and remove duplicates
# 6. Macs2: call peaks
# 7. IDR: extract reproducible peaks in replicates
# 8. Deeptools: create bigwig files

### ---------------
### Running command
### ---------------

# nohup bash ChipDNAv1.sh -c ChipDNAv1.conf &

### --------------
### Parse augments
### --------------

Help()
{
cat >&2 <<EOF
FUNCTION:
        Custome pipeline to perform ChIP-seq/CUT&RUN/CUT&Tag data analysis;
USAGE:
        $0 -c <configuration file> -h
        -c <*.conf file>      : configuration file that record all files needed by this pipeline;
        -h <help document>    : help document;
CONF FILE FORMAT:
        1. workding directory;
        2. full path of raw data;
        3. txt file that records sample metadata;
        4. data type: support PE and SE data;
        5. experimental design: data with input sample or not;
        6. target: support Histone and TFs;
        7. genomic sequence *.fa file;
        8. genome annotation *gtf file;
        9. adapter type: illumina/nextera;
        10. number of threads to use;
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
wd=`grep -v "^#" $cf | awk '{if(NR==1) print $0}'`;         echo -e "[1] working directory: ${wd}"
# 2. full path of raw data
rawdata=`grep -v "^#" $cf | awk '{if(NR==2) print $0}'`;    echo -e "[2] directory of raw fastq files: ${rawdata}"
# 3. txt file that records sample metadata
spanno=`grep -v "^#" $cf | awk '{if(NR==3) print $0}'`;     echo -e "[3] sample file: ${spanno}"
# 4. data type: support PE and SE data
dt=`grep -v "^#" $cf | awk '{if(NR==4) print $0}'`;         echo -e "[4] sample file: ${dt}"
# 5. experimental design: data with input sample or not
ed=`grep -v "^#" $cf | awk '{if(NR==5) print $0}'`;         echo -e "[5] data with input sample: ${ed}"
# 6. target: support Histone and TFs
tg=`grep -v "^#" $cf | awk '{if(NR==6) print $0}'`;         echo -e "[6] target: ${tg}"
# 7. genomic sequence *.fa file
fa=`grep -v "^#" $cf | awk '{if(NR==7) print $0}'`;         echo -e "[7] genome sequence file: ${fa}"
# 8. genome annotation *gtf file
gtf=`grep -v "^#" $cf | awk '{if(NR==8) print $0}'`;         echo -e "[8] genome annotation file: ${gtf}"
# 9. adapter type: illumina/nextera
adapter=`grep -v "^#" $cf | awk '{if(NR==9) print $0}'`;    echo -e "[9] adapter type: ${adapter}"
# 10. number of threads to use
nt=`grep -v "^#" $cf | awk '{if(NR==10) print $0}'`;        echo -e "[10] number of threads to use: ${nt}"
echo "> ---------------------------- <"

### --------------------
### Define the functions
### --------------------

# Create soft links for raw data
DataLink()
{
   indir=$1; spanno=$2; outdir=$3
   [ ! -d ${outdir} ] && mkdir -p ${outdir}
   cut -f 1 ${spanno} > ${wd}/metadata/sample_name.txt
   ls ${indir} | grep -f ${wd}/metadata/sample_name.txt - | while read file
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
   indir=$1; file=$2; outdir=$3; logdir=$4; qual=$5; len=$6; ada=$7; nt=$8
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
                     --${ada} --cores ${nt} > ${logdir}/${r1pre}_${r2pre}.log 2>&1
     done
   elif [ ${dt} == "SE" ]; then
     while read file
     do
         prefix=`echo ${file} | xargs basename -s ".fq.gz"`
         trim_galore ${indir}/${file} -o ${outdir} --quality ${qual} --max_n 4 \
                     --length ${len} --${ada} --cores ${nt} > ${logdir}/${prefix}.log 2>&1
     done < ${file}
   else
     echo "Unrecognized data type!"; exit 1
   fi
}
# Bowtie2 
# Build index
Bt2Index(){
   fa=$1; outdir=$2; logdir=$3; prefix=$4; nt=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   bowtie2-build --threads ${nt} ${fa} ${outdir}/${prefix} > ${logdir}/${prefix}.bt2.build.index.log 2>&1
}
# Mapping (customize the parameters by actual situation: Mismatch = 1 & no-mixed & no-discordant ==> by default)
Bt2Align(){
   indir=$1; file=$2; outdir=$3; logdir=$4; gx=$5; nt=$6
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ "${dt}" == "PE" ]; then
     for i in `seq 1 2 $(cat ${file} | wc -l)`
     do
         r1=$(awk -v row=${i} '(NR == row){print $0}' ${file})
         r1name=`echo ${r1} | xargs basename -s ".fq.gz"`
         r2=$(awk -v row=${i} '(NR == row+1){print $0}' ${file})
         r2name=`echo ${r2} | xargs basename -s ".fq.gz"`
         bowtie2 -t -q -N 1 -L 25 -p ${nt} --no-mixed --no-discordant -x ${gx} -1 ${indir}/${r1} -2 ${indir}/${r2} \
                 -S ${outdir}/${r1name}_${r2name}.sam > ${logdir}/${r1name}_${r2name}.log 2>&1
     done
   elif [ "${dt}" == "SE" ]; then
     cat ${file} | while read sam
     do
         name=`echo ${sam} | xargs basename -s ".fq.gz"`
         bowtie2 -t -q -N 1 -L 25 -p ${nt} -x ${gx} -U ${indir}/${sam} \
                 -S ${outdir}/${name}.sam > ${logdir}/${name}.log 2>&1
     done
   else
     echo "Error: Unrecognized data type!"; exit 1
   fi
}
# Sambamba
DedupSam()
{
   indir=$1; outdir=$2; np=$3; nt=$4
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}
   find ${indir} -name "*.sam" | xargs basename -s ".sam" | xargs -P ${np} -I{} sambamba view -t ${nt} -f bam -S -o ${outdir}/{}.bam ${indir}/{}.sam
   find ${outdir} -name "*.bam" | xargs basename -s ".bam" | xargs -P ${np} -I{} sambamba sort -o ${outdir}/{}_psort.bam ${outdir}/{}.bam
   find ${outdir} -name "*psort.bam" | xargs basename -s ".bam" | xargs -P ${np} -I{} sambamba markdup -t ${nt} -r ${outdir}/{}.bam ${outdir}/{}_dedup.bam
   find ${outdir} -type f ! -name '*psort*' -delete
}
# Macs2
Macs2HisPks()
{
   indir=$1; file=$2; outdir=$3; logdir=$4; mode=$5; pvalue=$6; spc=$7; ed=$8; dt=$9
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ "${ed}" == "T" ]; then
      awk '$2~"_IP$"' ${file} | cut -f 1,3 > ${wd}/metadata/call_pks.txt
      cat ${wd}/metadata/call_pks.txt | while read ip input
      do
          b1=`ls ${indir} | grep "${ip}" | grep "_dedup.bam$"`
          b2=`ls ${indir} | grep "${input}" | grep "_dedup.bam$"`
          name=$(basename -s ".bam" ${b1})
          if [ "${mode}" == "broad" ]; then
             if [ "${dt}" == "PE" ]; then
                macs2 callpeak -t ${indir}/${b1} -c ${indir}/${b2} -n ${name} -p ${pvalue} --broad --broad-cutoff ${pvalue} \
                               -g ${spc} -f BAMPE --outdir ${outdir} > ${logdir}/${name}.macs2.log 2>&1
             elif [ "${dt}" == "SE" ]; then
                macs2 callpeak -t ${indir}/${b1} -c ${indir}/${b2} -n ${name} -p ${pvalue} --broad --broad-cutoff ${pvalue} \
                               -g ${spc} --nomodel --shift 37 --extsize 73 --outdir ${outdir} > ${logdir}/${name}.macs2.log 2>&1
             else 
                echo "Error: Unrecognized data type!"; exit 1
             fi
          elif [ "${mode}" == "narrow" ]; then
             if [ "${dt}" == "PE" ]; then
                macs2 callpeak -t ${indir}/${b1} -c ${indir}/${b2} -n ${name} -p ${pvalue} \
                               -g ${spc} -f BAMPE --outdir ${outdir} > ${logdir}/${name}.macs2.log 2>&1
             elif [ "${dt}" == "SE" ]; then
                macs2 callpeak -t ${indir}/${b1} -c ${indir}/${b2} -n ${name} -p ${pvalue} \
                               -g ${spc} --nomodel --shift 37 --extsize 73 --outdir ${outdir} > ${logdir}/${name}.macs2.log 2>&1
             else
                echo "Error: Unrecognized data type!"; exit 1
             fi
          else
             echo "Error: Unrecognized peak type!"; exit 1
          fi
      done
   elif [ "${ed}" == "F" ]; then
      awk '$2~"_IP$"||$2~"_IgG$"' ${file} | cut -f 1 > ${wd}/metadata/call_pks.txt
      cat ${wd}/metadata/call_pks.txt | while read ip
      do  
          b1=`ls ${indir} | grep "${ip}" | grep "_dedup.bam$"`
          name=$(basename -s ".bam" ${b1})
          if [ "${mode}" == "broad" ]; then
             if [ "${dt}" == "PE" ]; then
                macs2 callpeak -t ${indir}/${b1} -n ${name} -p ${pvalue} --nolambda --broad --broad-cutoff ${pvalue} \
                               -g ${spc} -f BAMPE --outdir ${outdir} > ${logdir}/${name}.macs2.log 2>&1
             elif [ "${dt}" == "SE" ]; then
                macs2 callpeak -t ${indir}/${b1} -n ${name} -p ${pvalue} --nolambda --broad --broad-cutoff ${pvalue} \
                               -g ${spc} --nomodel --shift 37 --extsize 73 --outdir ${outdir} > ${logdir}/${name}.macs2.log 2>&1
             else 
                echo "Error: Unrecognized data type!"; exit 1
             fi
          elif [ "${mode}" == "narrow" ]; then
             if [ "${dt}" == "PE" ]; then
                macs2 callpeak -t ${indir}/${b1} -n ${name} -p ${pvalue} --nolambda \
                               -g ${spc} -f BAMPE --outdir ${outdir} > ${logdir}/${name}.macs2.log 2>&1
             elif [ "${dt}" == "SE" ]; then
                macs2 callpeak -t ${indir}/${b1} -n ${name} -p ${pvalue} --nolambda \
                               -g ${spc} --nomodel --shift 37 --extsize 73 --outdir ${outdir} > ${logdir}/${name}.macs2.log 2>&1
             else
                echo "Error: Unrecognized data type!"; exit 1
             fi
          else
             echo "Error: Unrecognized peak type!"; exit 1
          fi
      done
   else
      echo "Error: Unrecognized experimental design type!"; exit 1
   fi
}
Macs2TFsPks()
{
   indir=$1; file=$2; outdir=$3; logdir=$4; pvalue=$5; spc=$6; ed=$7; dt=$8
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ "${ed}" == "T" ]; then
      awk '$2~"_IP$"' ${file} | cut -f 1,3 > ${wd}/metadata/call_pks.txt
      cat ${wd}/metadata/call_pks.txt | while read ip input
      do
          b1=`ls ${indir} | grep "${ip}" | grep "_dedup.bam$"`
          b2=`ls ${indir} | grep "${input}" | grep "_dedup.bam$"`
          name=$(basename -s ".bam" ${b1})
          if [ "${dt}" == "PE" ]; then
             macs2 callpeak -t ${indir}/${b1} -c ${indir}/${b2} -n ${name} -p ${pvalue} \
                            -g ${spc} -f BAMPE --outdir ${outdir} > ${logdir}/${name}.macs2.log 2>&1
          elif [ "${dt}" == "SE" ]; then
             macs2 callpeak -t ${indir}/${b1} -c ${indir}/${b2} -n ${name} -p ${pvalue} \
                            -g ${spc} --nomodel --shift 100 --extsize 200 --outdir ${outdir} > ${logdir}/${name}.macs2.log 2>&1
          else
             echo "Error: Unrecognized data type!"; exit 1
          fi
      done
   elif [ "${ed}" == "F" ]; then
      awk '$2~"_IP$"||$2~"_IgG$"' ${file} | cut -f 1 > ${wd}/metadata/call_pks.txt
      cat ${wd}/metadata/call_pks.txt | while read ip
      do
          b1=`ls ${indir} | grep "${ip}" | grep "_dedup.bam$"`
          name=$(basename -s ".bam" ${b1})
          if [ "${dt}" == "PE" ]; then
             macs2 callpeak -t ${indir}/${b1} -n ${name} -p ${pvalue} --nolambda \
                            -g ${spc} -f BAMPE --outdir ${outdir} > ${logdir}/${name}.macs2.log 2>&1
          elif [ "${dt}" == "SE" ]; then
             macs2 callpeak -t ${indir}/${b1} -n ${name} -p ${pvalue} --nolambda \
                            -g ${spc} --nomodel --shift 100 --extsize 200 --outdir ${outdir} > ${logdir}/${name}.macs2.log 2>&1
          else
             echo "Error: Unrecognized data type!"; exit 1
          fi
      done
   else
      echo "Error: Unrecognized experimental design type!"; exit 1
   fi
}
# idr
IDR(){
   indir=$1; file=$2; outdir=$3; logdir=$4; mode=$5
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   pk1=$(awk 'NR == 1' ${file}); pk1_pre="${pk1%%_psort_dedup*}"; pk1_post="${pk1##*.}"
   pk2=$(awk 'NR == 2' ${file}); pk2_pre="${pk2%%_psort_dedup*}"; pk2_post="${pk2##*.}"
   sort -k8,8nr ${indir}/${pk1} > ${outdir}/${pk1_pre}_sorted.${pk1_post}
   sort -k8,8nr ${indir}/${pk2} > ${outdir}/${pk2_pre}_sorted.${pk2_post}
   idr --samples ${outdir}/${pk1_pre}_sorted.${pk1_post} ${outdir}/${pk2_pre}_sorted.${pk2_post} \
       --input-file-type ${mode} --rank p.value --output-file ${outdir}/${pk1_pre}_${pk2_pre}.${mode} \
       --plot --log-output-file ${logdir}/${pk1_pre}_${pk2_pre}.log
}
# Deeptools
DeepCoverage(){
   indir=$1; file=$2; outdir=$3; logdir=$4; gs=$5; nt=$6
   [ ! -d "${outdir}" ] && mkdir -p ${outdir}; [ ! -d "${logdir}" ] && mkdir -p ${logdir}
   if [ ${dt} == "PE" ]; then
      cat ${file} | while read bam
      do
          prefix=$(basename -s ".bam" ${bam})
          bamCoverage --bam ${indir}/${bam} --outFileName ${outdir}/${prefix}.bw --outFileFormat bigwig \
                      --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize ${gs} \
                      --extendReads --ignoreDuplicates --numberOfProcessors ${nt} > ${logdir}/${prefix}_coverage.log 2>&1
      done
   elif [ ${dt} == "SE" ]; then
      cat ${file} | while read bam
      do
          prefix=$(basename -s ".bam" ${bam})
          bamCoverage --bam ${indir}/${bam} --outFileName ${outdir}/${prefix}.bw --outFileFormat bigwig \
                      --binSize 10 --normalizeUsing RPKM --effectiveGenomeSize ${gs} \
                      --ignoreDuplicates --numberOfProcessors ${nt} > ${logdir}/${prefix}_coverage.log 2>&1
      done
   else
      echo "Error: Unrecognized data type!"; exit 1
   fi
}

### ---
### Run
### ---

echo "1st step. Build working shop"
echo "Begin at $(date)"
mkdir logs results rawdata metadata scripts
DataLink ${rawdata} ${spanno} ${wd}/rawdata
[ ! -e "${wd}/scripts/faCount" ] && wget -O ${wd}/scripts/faCount http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faCount
chmod +x ${wd}/scripts/faCount
${wd}/scripts/faCount ${fa} > ${wd}/metadata/genome_size.txt
egs=`awk 'NR>1' ${wd}/metadata/genome_size.txt | grep -v "total" | awk 'BEGIN{total=0}{total+=($2-$7)}END{print total}'`
echo "The effective genome size is ${egs}"
echo "Finish at $(date)"



echo "2nd step. Check the quality of raw data"
echo "Begin at $(date)"
QcReads ${wd}/rawdata ${wd}/results/fastqc/rawdata ${wd}/logs/fastqc/rawdata
ParseQc ${wd}/results/fastqc/rawdata ${wd}/results/multiqc/rawdata ${wd}/logs/multiqc/rawdata
ls ${wd}/rawdata | grep "fq.gz$" > ${wd}/metadata/raw_fq_${dt}.txt
echo "Finish at $(date)"



echo "3rd step. Filter and trim the reads"
echo "Begin at $(date)"
quality=20; length=30
TrimReads ${wd}/rawdata ${wd}/metadata/raw_fq_${dt}.txt ${wd}/results/trimgalore ${wd}/logs/trimgalore ${quality} ${length} ${adapter} ${nt}
if [ "${dt}" == "PE" ]; then
   ls ${wd}/results/trimgalore | grep "fq.gz$" | grep "val" > ${wd}/metadata/trimgalore_fq_${dt}.txt
else
   ls ${wd}/results/trimgalore | grep "fq.gz$" | grep "trimmed" > ${wd}/metadata/trimgalore_fq_${dt}.txt
fi
echo "Finish at $(date)"



echo "4th step. Check the quality of raw data again"
echo "Begin at $(date)"
QcReads ${wd}/results/trimgalore ${wd}/results/fastqc/trimgalore ${wd}/logs/fastqc/trimgalore
ParseQc ${wd}/results/fastqc/trimgalore ${wd}/results/multiqc/trimgalore ${wd}/logs/multiqc/trimgalore
echo "Finish at $(date)"



echo "5th step. Align read to genome"
echo "Begin at $(date)"
prefix=$(basename -s ".fa" ${fa})
Bt2Index ${fa} ${wd}/results/bowtie2/index ${wd}/logs/bowtie2/index ${prefix} ${nt}
Bt2Align ${wd}/results/trimgalore ${wd}/metadata/trimgalore_fq_${dt}.txt ${wd}/results/bowtie2/align ${wd}/logs/bowtie2/align ${wd}/results/bowtie2/index/${prefix} ${nt}
ParseQc ${wd}/logs/bowtie2/align ${wd}/results/multiqc/bowtie2 ${wd}/logs/multiqc/bowtie2
if [ "${dt}" == "PE" ]; then
   ls ${wd}/results/bowtie2/align | grep "val" | grep "sam$" > ${wd}/metadata/sam_${dt}.txt
else
   ls ${wd}/results/bowtie2/align | grep "trimmed" | grep "sam$" > ${wd}/metadata/sam_${dt}.txt
fi
echo "Finish at $(date)"



echo "6th step. Process the sam files"
echo "Begin at $(date)"
if [ $(wc -l ${wd}/metadata/sam_${dt}.txt | awk '{print $1}') -ge 3 ]; then
   np=3
else
   np=1
fi
DedupSam ${wd}/results/bowtie2/align ${wd}/results/sambamba ${np} ${nt}
if [ "${dt}" == "PE" ]; then
   ls ${wd}/results/sambamba | grep "val" | grep "dedup.bam$" > ${wd}/metadata/psort_dedup_bam_${dt}.txt
else
   ls ${wd}/results/sambamba | grep "trimmed" | grep "dedup.bam$" > ${wd}/metadata/psort_dedup_bam_${dt}.txt
fi
echo "Finish at $(date)"



echo "7th step. Call peaks"
echo "Begin at $(date)"
if [ "${tg}" == "Histone" ]; then
   for pt in broad narrow
   do
      for p in 0.05 0.01 0.005
      do
         echo "Calling ${pt} peaks for histone modification data"
         Macs2HisPks ${wd}/results/sambamba ${wd}/sample_anno.txt ${wd}/results/macs2/his_${pt}Peak_p${p} \
                     ${wd}/logs/macs2/his_${pt}Peak_p${p} ${pt} ${p} ${egs} ${ed} ${dt}
      done
   done
elif [ "${tg}" == "TFs" ]; then
   echo "Calling narrow peaks for transcription factors data"
   for p in 0.05 0.01 0.005
   do
      Macs2TFsPks ${wd}/results/sambamba ${wd}/sample_anno.txt ${wd}/results/macs2/tfs_narrowPeak_p${p} \
                  ${wd}/logs/macs2/tfs_narrowPeak_p${p} ${p} ${egs} ${ed} ${dt}
   done
else
   echo "Error! Unrecognized data target type!"; exit 1
fi
echo "Finish at $(date)"



echo "8th step. IDR"
echo "Begin at $(date)"
awk '$2~"_IP$"' ${spanno} | cut -f 2 | sort -u | while read group
do
   awk -v g=${group} '$2==g' ${spanno} | cut -f 1 > ${wd}/metadata/${group}_sp_prefix.txt
   n=`wc -l ${wd}/metadata/${group}_sp_prefix.txt | awk '{print $1}'`
   if [ "${n}" -le 1 ]; then
      echo "Only one sample for group ${group}..."
   elif [ "${n}" -gt 1 ]; then
      ls ${wd}/results/macs2 | while read dir
      do
         mode=$(echo "${dir}" | cut -d "_" -f 2)
         pvalue=$(echo "${dir}" | cut -d "_" -f 3)
         echo "Run IDR for group ${group} in ${mode} mode under ${pvalue}"
         ls ${wd}/results/macs2/${dir} | grep "Peak$" | grep -f ${wd}/metadata/${group}_sp_prefix.txt - > ${wd}/metadata/${group}_${mode}_idr.txt
         IDR ${wd}/results/macs2/${dir} ${wd}/metadata/${group}_${mode}_idr.txt ${wd}/results/idr/${dir} ${wd}/logs/idr/${dir} ${mode}
      done
   fi
done
echo "Finish at $(date)"



echo "9th step. Creating bigwig files"
echo "Begin at $(date)"
DeepCoverage ${wd}/results/sambamba ${wd}/metadata/psort_dedup_bam_${dt}.txt ${wd}/results/deeptools/coverage \
             ${wd}/logs/deeptools/coverage ${egs} ${nt}
echo "Finish at $(date)"
