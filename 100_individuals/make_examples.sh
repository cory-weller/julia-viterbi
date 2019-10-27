# !/bin/bash

echo "Retrieving the .vcf files"
wget -L https://virginia.box.com/shared/static/d9943yx9let83lwm7aipn71m83ylztba.gz -O 32_F5_1-100.tar.gz

echo "Extracting"
tar -xvzf 32_F5_1-100.tar.gz

python scripts/reformat_hap_files.py *.vcf

echo "Retrieving the haplotypes.polarized.vcf file"
wget -L https://virginia.box.com/shared/static/4r7uqcrk00dlinsc6134yte22vnjvy53.gz -O haplotypes.polarized.vcf.gz

echo "Extracting"
gunzip haplotypes.polarized.vcf.gz

echo "Reformatting haplotypes.polarized.vcf to fit our notation standards"
sed -i 's/\t0/\t./g' haplotypes.polarized.vcf
echo "1/3 done"
sed -i 's/\t1/\t0/g' haplotypes.polarized.vcf
echo "2/3 done"
sed -i 's/\t-1/\t1/g' haplotypes.polarized.vcf
echo "Done"

python scripts/reformat_hap_files.py *.wide.haps

echo "Downsampling"
scripts/downsample.sh

#Now I think we just need to make the strips and move them in the right place...
