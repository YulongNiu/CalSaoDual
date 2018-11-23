date

BEDTOOLS_PATH=/home/Yulong/Biotools/bedtools2/bin

INDEX_PATH=/home/Yulong/Biotools/RefData/sao

cd ${INDEX_PATH}

${BEDTOOLS_PATH}/bedtools getfasta -fi NC_007795.fna -bed NC_007795.gff -name -s -fullHeader -fo NC_007795_cdna.fa

date
