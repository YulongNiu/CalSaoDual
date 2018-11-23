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

cd ${INDEX_PATH}

## pseudoquant
${KALLISTO_PATH}/kallisto quant \
                -i cal_sao_idx \
                -t 8 \
                -o testpseudo \
                -l \
                cal_sao_R1.fq.gz cal_sao_R2.fq.gz


# ## quantification
# cd ${CLEAN_PATH}


# fmember=('4h_sm_1' '4h_sm_2' '4h_sm_3' '4h_srta_1' '4h_srta_2' '4h_srta_3' '24h_sm_1' '24h_sm_2' '24h_sm_3' '24h_srta_1' '24h_srta_2' '24h_srta_3')

# for i in "${fmember[@]}"
# do
#     echo "Quantify ${i}"
#     ${KALLISTO_PATH}/kallisto quant \
#                     -t 8 \
#                     -i ${INDEX_PATH}/smu_kallisto_idx \
#                     -o ${RES_PATH}/${i} \
#                     ${i}_R1.fq.gz ${i}_R2.fq.gz
# done

date
