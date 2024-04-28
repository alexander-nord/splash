
OUT_DIR_NAME='genomes'

if [[ -d $OUT_DIR_NAME ]]
then
    rm -rf $OUT_DIR_NAME
fi

mkdir $OUT_DIR_NAME
cd $OUT_DIR_NAME

wget 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz' -O hg38.fa.gz
gunzip hg38.fa.gz

wget 'https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz' -O mm39.fa.gz
gunzip mm39.fa.gz

wget 'https://hgdownload.soe.ucsc.edu/goldenPath/rn6/bigZips/rn6.fa.gz' -O rn6.fa.gz
gunzip rn6.fa.gz

echo "homo_sapiens      hg38.fa"  > "species-to-genome"
echo "mus_musculus      mm39.fa" >> "species-to-genome"
echo "rattus_norvegicus rn6.fa"  >> "species-to-genome"
echo "human             hg38.fa" >> "species-to-genome"
echo "mouse             mm39.fa" >> "species-to-genome"
echo "rat               rn6.fa"  >> "species-to-genome"

