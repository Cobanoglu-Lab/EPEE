# # paramscan test
# python paramscan.py -a ../../data/rnaseq/immune/CD4_Naive.txt.gz -b ../../data/rnaseq/immune/CD4_Th2.txt.gz -na ../../data/network/cd4+_t_cells.txt.gz -nb ../../data/network/cd4+_t_cells.txt.gz -o ../../output/paramscan/Th2/ -eval Th2 -r 2 -i 100 -prefix Th2 -nodenumbers 4 -multiprocess


# # AML
# python paramscan.py -a ../../data/microarray/aml/AML_normal.txt.gz -b ../../data/microarray/aml/AML_cancer.txt.gz -na ../../data/network/11_myeloid_leukocytes.txt.gz -nb ../../data/network/15_myeloid_leukemia.txt.gz -o ../../output/paramscan/AML/ -eval AML -prefix AML -nodenumbers 4 -multiprocess -run

# Th1
python paramscan.py -a ../../data/rnaseq/immune/CD4_Naive.txt.gz -b ../../data/rnaseq/immune/CD4_Th1.txt.gz -na ../../data/network/cd4+_t_cells.txt.gz -nb ../../data/network/cd4+_t_cells.txt.gz -o ../../output/paramscan/Th1/ -eval Th1 -prefix Th1 -nodenumbers 4 -multiprocess -run

# Th2
python paramscan.py -a ../../data/rnaseq/immune/CD4_Naive.txt.gz -b ../../data/rnaseq/immune/CD4_Th2.txt.gz -na ../../data/network/cd4+_t_cells.txt.gz -nb ../../data/network/cd4+_t_cells.txt.gz -o ../../output/paramscan/Th2/ -eval Th2 -prefix Th2 -nodenumbers 4 -multiprocess -run

# Th17
python paramscan.py -a ../../data/rnaseq/immune/CD4_Naive.txt.gz -b ../../data/rnaseq/immune/CD4_Th17.txt.gz -na ../../data/network/cd4+_t_cells.txt.gz -nb ../../data/network/cd4+_t_cells.txt.gz -o ../../output/paramscan/Th2/ -eval Th17 -prefix Th17 -nodenumbers 4 -multiprocess -run

# Bmem
python paramscan.py -a ../../data/rnaseq/immune/B_Naive.txt.gz -b ../../data/rnaseq/immune/B_Memory.txt.gz -na ../../data/network/cd8+_t_cells.txt.gz -nb ../../data/network/cd8+_t_cells.txt.gz -o ../../output/paramscan/Bmem/ -eval Bmem -prefix Bmem -nodenumbers 4 -multiprocess -run

# COAD
python paramscan.py -a ../../data/rnaseq/tcga/TCGA_COAD_SolidTissueNormal_FPKM_UQ.txt.gz -b ../../data/rnaseq/tcga/TCGA_COAD_PrimaryTumor_FPKM_UQ.txt.gz -na ../../data/network/20_gastrointestinal_system.txt.gz -nb ../../data/network/20_gastrointestinal_system.txt.gz -o ../../output/paramscan/COAD/ -eval COAD -prefix COAD -nodenumbers 4 -multiprocess -run
