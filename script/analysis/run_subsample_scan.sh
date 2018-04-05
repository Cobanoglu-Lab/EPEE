# # Th2
python subsample_scan.py -a ../../data/rnaseq/immune/CD4_Naive.txt.gz -b ../../data/rnaseq/immune/CD4_Th2.txt.gz -o ../../data/rnaseq/ -prefix Th2
# sbatch ../sruns/20180223_Th2_subsamples.sh

# COAD
# python subsample_scan.py -a ../../data/rnaseq/tcga/TCGA_COAD_SolidTissueNormal_FPKM_UQ.txt.gz -b ../../data/rnaseq/tcga/TCGA_COAD_PrimaryTumor_FPKM_UQ.txt.gz -o ../../data/rnaseq/ -prefix COAD
# sbatch ../sruns/20180223_COAD_subsamples.sh
