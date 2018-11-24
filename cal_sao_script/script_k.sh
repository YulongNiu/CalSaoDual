date

KALLISTO_PATH=/home/Yulong/Biotools/kallisto_linux-v0.44.0

PROJECT_PATH=/extDisk2/cal_sao

INDEX_PATH=${PROJECT_PATH}/kallisto_results
CLEAN_PATH=${PROJECT_PATH}/cleandata

## merge fa
# cat ${INDEX_PATH}/cal_cdna_name.fa ${INDEX_PATH}/sao_cdna_name.fa > ${INDEX_PATH}/cal_sao_cdna_name.fa

# ## build index
# ${KALLISTO_PATH}/kallisto index \
#                 -i ${INDEX_PATH}/cal_sao_idx \
#                 ${INDEX_PATH}/cal_sao_cdna_name.fa

# ${KALLISTO_PATH}/kallisto index \
#                 -i ${INDEX_PATH}/cal_idx \
#                 ${INDEX_PATH}/cal_cdna_name.fa

# ${KALLISTO_PATH}/kallisto index \
#                 -i ${INDEX_PATH}/sao_idx \
#                 ${INDEX_PATH}/sao_cdna_name.fa


cd ${CLEAN_PATH}

# ## pseudoquant
# ${KALLISTO_PATH}/kallisto quant \
#                 -i ${INDEX_PATH}/cal_sao_idx \
#                 -t 8 \
#                 -o testpseudo \
#                 Sa_CAF2_1_1_R1.fq.gz Sa_CAF2_1_1_R2.fq.gz


## quantification
## cal
# callabel=('CAF2_1_1' 'CAF2_1_2' 'CAF2_1_3')
# for i in "${callabel[@]}"
# do
#     echo "Quantify ${i}"
#     ${KALLISTO_PATH}/kallisto quant \
#                     -t 8 \
#                     -i ${INDEX_PATH}/cal_idx \
#                     -o ${INDEX_PATH}/${i} \
#                     ${i}_R1.fq.gz ${i}_R2.fq.gz
# done

# ## sao
# saolabel=('Sa_1' 'Sa_2' 'Sa_3')
# for i in "${saolabel[@]}"
# do
#     echo "Quantify ${i}"
#     ${KALLISTO_PATH}/kallisto quant \
#                     -t 8 \
#                     -i ${INDEX_PATH}/sao_idx \
#                     -o ${INDEX_PATH}/${i} \
#                     ${i}_R1.fq.gz ${i}_R2.fq.gz
# done

## cal_sao
calsaolabel=('Sa_CAF2_1_1' 'Sa_CAF2_1_2' 'Sa_CAF2_1_3')
for i in "${calsaolabel[@]}"
do
    echo "Quantify ${i}"
    ${KALLISTO_PATH}/kallisto quant \
                    -t 8 \
                    -i ${INDEX_PATH}/cal_sao_idx \
                    -o ${INDEX_PATH}/${i} \
                    ${i}_R1.fq.gz ${i}_R2.fq.gz

    ## pseudoquant
    echo "Pseudo ${i}"
    ${KALLISTO_PATH}/kallisto pseudo \
                        -i ${INDEX_PATH}/cal_sao_idx \
                        -o ${INDEX_PATH}/${i}_pseudo \
                        ${i}_R1.fq.gz ${i}_R2.fq.gz
done

date
