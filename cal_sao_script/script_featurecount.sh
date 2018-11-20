date

FEATURECOUNT_PATH=/home/Yulong/Biotools/subread-1.6.2-Linux-x86_64/bin

CALGFF_PATH=/home/Yulong/Biotools/RefData/cal
SAOGFF_PATH=/home/Yulong/Biotools/RefData/sao

PROJECT_PATH=/extDisk2/cal_sao
FILTER_PATH=${PROJECT_PATH}/filter_alignments
SAVE_PATH=${PROJECT_PATH}/figures_tables

CORENUM=8

cd ${FILTER_PATH}

## cal counts
${FEATURECOUNT_PATH}/featureCounts \
                    -p \
                    -T ${CORENUM} \
                    -a ${CALGFF_PATH}/C_albicans_SC5314_A22_current_features.gff \
                    -M \
                    -t gene \
                    -g ID \
                    -o ${SAVE_PATH}/cal_counts.txt \
                    `ls | grep 'ca_sort.bam$'`

## sao counts
${FEATURECOUNT_PATH}/featureCounts \
                    -p \
                    -T ${CORENUM} \
                    -a ${SAOGFF_PATH}/Sa.gff \
                    -M \
                    -t gene \
                    -g ID \
                    -o ${SAVE_PATH}/sao_counts.txt \
                    `ls | grep 'sa_sort.bam$'`

date
