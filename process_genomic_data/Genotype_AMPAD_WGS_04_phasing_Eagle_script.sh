wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.1.tar.gz
tar xzvf Eagle_v2.4.1.tar.gz 

mkdir input
cd input
gsutil cp gs://shinya_test/Genotype/AMPAD_WGS/01_qc/output/pca_outliers.txt ./
  
PLINK_HOME=/data/shinya/temp/
BFILE=../../temp/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08.recalibrated_variants_hg38_dn.QCed1
$PLINK_HOME/plink --bfile $BFILE  --remove pca_outliers.txt  --keep-allele-order --split-x hg38 --make-bed --out plink
for i in `seq 1 23`;
do
#$PLINK_HOME/plink --bfile plink  --remove pca_outliers.txt --keep-allele-order --real-ref-alleles --chr $i --maf 0.0004 --hwe 1e-50 --geno 0.1  --recode vcf  --out chr$i
bgzip -c chr$i.vcf >  chr$i.vcf.gz
done    

for i in $(seq 1 23); do
echo "../Eagle_v2.4.1/eagle --vcf=chr$i.vcf --geneticMapFile=../Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz --numThreads 8 --vcfOutFormat=z --outPrefix=../output/chr$i.phased > ../output/chr$i.phased.log" >> myCommands.txt
done

mkdir ../output
less myCommands.txt | parallel -j8 {}


# combine
/home/shinya/Resource/GENETICS_Resource/softwares/bcftools-1.10.2/bcftools concat -Oz /data/shinya/temp_eagle/output/*.vcf.gz > /data/shinya/temp_eagle/output/chrAll.phased.vcf.gz
/home/shinya/Resource/GENETICS_Resource/softwares/bcftools-1.10.2/bcftools index /data/shinya/temp_eagle/output/chrAll.phased.vcf.gz --threads 8
~/Resource/NGS_Resource/softwares/htslib-1.9/tabix -p vcf /data/shinya/temp_eagle/output/chrAll.phased.vcf.gz

# plink
PLINK_HOME=/data/shinya/temp/
$PLINK_HOME/plink --vcf /data/shinya/temp_eagle/output/chrAll.phased.vcf.gz --keep-allele-order --make-bed --out /data/shinya/temp_eagle/output/chrAll.phased

