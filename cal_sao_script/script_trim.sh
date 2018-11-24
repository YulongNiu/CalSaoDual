date

BASE_PATH=/extDisk2/cal_sao
RAW_PATH=${BASE_PATH}/rawdata/baidu_raw
CLEAN_PATH=${BASE_PATH}/cleandata

FASTP_PATH=/home/Yulong/Biotools/fastp

CORENUM=8

rawlabel=('CAF2_1_1' 'CAF2_1_2' 'CAF2_1_3' 'Sa_1' 'Sa_2' 'Sa_3' 'Sa_CAF2_1_1' 'Sa_CAF2_1_2' 'Sa_CAF2_1_3')
# cleanlabel=('cal_1' 'cal_2' 'cal_3' 'sao_1' 'sao_2' 'sao_3' 'cal_sao_1' 'cal_sao_2' 'cal_sao_3')


cd ${RAW_PATH}

for i in "${rawlabel[@]}"
do
    eachraw=($(ls -d -1 ${RAW_PATH}/${i}/*))
    echo "Trimming ${eachraw[0]} and ${eachraw[1]}."

    ${FASTP_PATH}/fastp -w ${CORENUM} \
                 -z 6 \
                 -g -p -c \
                 -h ${i}.html -j ${i}.json \
                 -i ${eachraw[0]} -I ${eachraw[1]} \
                 -o ${CLEAN_PATH}/${i}_R1.fq.gz -O ${CLEAN_PATH}/${i}_R2.fq.gz

done

date
