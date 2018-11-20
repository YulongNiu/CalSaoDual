date

SAMTOOLS_PATH=/home/Yulong/Biotools/samtools_1.8/bin/

FILTER_PATH=/extDisk2/cal_sao/filter_alignments/

CORENUM=8

cd ${FILTER_PATH}

bamfiles=`ls | grep .bam`

for i in ${bamfiles}
do
    echo sort and index ${i}

    sortname=${i%.*}

    ## sort
    ${SAMTOOLS_PATH}/samtools sort \
                    -@ ${CORENUM} \
                    -o ${sortname}_sort.bam \
                    ${i}

    ## index
    ${SAMTOOLS_PATH}/samtools index \
                    ${sortname}_sort.bam

done

date
