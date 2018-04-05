# Th1
python run_epee_v0.1.4.3.py -a ../data/rnaseq/immune/CD4_Naive.txt.gz -b ../data/rnaseq/immune/CD4_Th1.txt.gz -na ../data/network/cd4+_t_cells.txt.gz -nb ../data/network/cd4+_t_cells.txt.gz -o ../output/ -prefix Th1

# Th2
python run_epee_v0.1.4.3.py -a ../data/rnaseq/immune/CD4_Naive.txt.gz -b ../data/rnaseq/immune/CD4_Th2.txt.gz -na ../data/network/cd4+_t_cells.txt.gz -nb ../data/network/cd4+_t_cells.txt.gz -o ../output/ -prefix Th2

# Th17
python run_epee_v0.1.4.3.py -a ../data/rnaseq/immune/CD4_Naive.txt.gz -b ../data/rnaseq/immune/CD4_Th17.txt.gz -na ../data/network/cd4+_t_cells.txt.gz -nb ../data/network/cd4+_t_cells.txt.gz -o ../output/ -prefix Th17

# Bmem
python run_epee_v0.1.4.3.py -a ../data/rnaseq/immune/B_Naive.txt.gz -b ../data/rnaseq/immune/B_Memory.txt.gz -na ../data/network/cd8+_t_cells.txt.gz -nb ../data/network/cd8+_t_cells.txt.gz -o ../output/ -prefix Bmem

# TCGA COAD
python run_epee_v0.1.4.3.py -a ../data/rnaseq/tcga/TCGA_COAD_SolidTissueNormal_FPKM_UQ.txt.gz -b ../data/rnaseq/tcga/TCGA_COAD_PrimaryTumor_FPKM_UQ.txt.gz -na ../data/network/20_gastrointestinal_system.txt.gz -nb ../data/network/20_gastrointestinal_system.txt.gz -o ../output/ -prefix COAD

# AML
python run_epee_v0.1.4.3.py -a ../data/microarray/aml/AML_normal.txt.gz -b ../data/microarray/aml/AML_cancer.txt.gz -na ../data/network/11_myeloid_leukocytes.txt.gz -nb ../data/network/15_myeloid_leukemia.txt.gz -o ../output/ -prefix AML

# Running different linear models to evaluate the performance of GCFL (Graph
# constrained fused lasso). By default EPEE runs GCFL model. No penalty model
# is the model that does not have either lasso or GC (Graph constrained). GC
# model is just with the graph constrained regularizer without lasso
# regularization.

# Th2 EPEE
python ../run_epee_v0.1.4.3.py -a ../../data/rnaseq/immune/CD4_Naive.txt.gz -b ../../data/rnaseq/immune/CD4_Th2.txt.gz -na ../../data/network/cd4+_t_cells.txt.gz -nb ../../data/network/cd4+_t_cells.txt.gz -o ../../output/Th2/ -eval Th2 -prefix Th2_epee -multiprocess

# Th2 EPEE no penalty
python ../run_epee_v0.1.4.3.py -a ../../data/rnaseq/immune/CD4_Naive.txt.gz -b ../../data/rnaseq/immune/CD4_Th2.txt.gz -na ../../data/network/cd4+_t_cells.txt.gz -nb ../../data/network/cd4+_t_cells.txt.gz -o ../../output/Th2/ -eval Th2 -prefix Th2_no_penalty -model no-penalty -multiprocess

# Th2 EPEE GC
python ../run_epee_v0.1.4.3.py -a ../../data/rnaseq/immune/CD4_Naive.txt.gz -b ../../data/rnaseq/immune/CD4_Th2.txt.gz -na ../../data/network/cd4+_t_cells.txt.gz -nb ../../data/network/cd4+_t_cells.txt.gz -o ../../output/Th2/ -eval Th2 -prefix Th2_GC -reg1 0 -reg2 0.01 -multiprocess

# Th2 EPEE lasso
python ../run_epee_v0.1.4.3.py -a ../../data/rnaseq/immune/CD4_Naive.txt.gz -b ../../data/rnaseq/immune/CD4_Th2.txt.gz -na ../../data/network/cd4+_t_cells.txt.gz -nb ../../data/network/cd4+_t_cells.txt.gz -o ../../output/Th2/ -eval Th2 -prefix Th2_lasso -model epee-l -multiprocess
