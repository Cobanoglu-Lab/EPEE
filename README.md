# EPEE

Effectors and Perturbation Estimation Engine (EPEE) a sparse linear model with graph constrained lasso regularization for differential analysis of RNA-seq data. The inputs are transcriptomic data for the two conditions under comparison, and context-specific TF-gene networks. If transcriptomic data is sequencing based, then data needs to be normalized to either TPM/FPKM/RPKM. EPEE is implemented in Python, using TensorFlow.

### Inputs

- conditionA.txt and conditionB.txt

EPEE requires expression data matrix for the two conditions. Input is a tab delimited file in which columns are the samples and rows are the genes. First column of the file needs to be gene names.
Please do not log normalize the dataset before running EPEE. EPEE log transforms the data to log(TPM/FPKM/RPKM + 1).

- networkA and networkB

EPEE requires context specific networks. Currently EPEE supports 426 context-specific networks published by Marbach et al. Nature Methods 2016.


### Setup

1. Install [Anaconda](https://www.anaconda.com/download)

2. Download the [networks](http://regulatorycircuits.org/download.html) and example [data](https://github.com/Cobanoglu-Lab/EPEE/tree/master/test/data) to run EPEE.

3. Clone the git repository and set up the conda environment to run EPEE. We provided the environment files in `env` directory. If your machine has GPU card then we recommend that you use `epee_GPU.txt` file to create environment, otherwise create environment using `epee_CPU.txt`. To create conda environment use following command
```
conda env create -f epee_CPU.txt -n epee
```
Activate the new environment: Windows: `activate epee`, macOS and Linux: `source activate epee`

4. View the available human networks, and determine the network appropriate for your context.

5. Usage to run EPEE

```
python run_epee.py -a <conditionA.txt>
                   -b <conditionB.txt>
                   -na <networkA.txt>
                   -nb <networkB.txt>
                   -o <output_directory>
```

### Example

##### CD4 Naive vs Th2 differential analysis

```
python run_epee_v0.1.4.3.py -a ../data/rnaseq/immune/CD4_Naive.txt.gz
                            -b ../data/rnaseq/immune/CD4_Th2.txt.gz
                            -na ../data/network/cd4+_t_cells.txt.gz
                            -nb ../data/network/cd4+_t_cells.txt.gz
                            -o /path/to/output_directory/
                            -prefix Th2
```

##### Normal Colon vs Colorectal Adenocarcinoma (COAD) differential analysis

```
python run_epee_v0.1.4.3.py -a ../data/rnaseq/tcga/TCGA_COAD_SolidTissueNormal_FPKM_UQ.txt
                            -b ../data/rnaseq/tcga/TCGA_COAD_PrimaryTumor_FPKM_UQ.txt
                            -na ../data/network/20_gastrointestinal_system.txt.gz
                            -nb ../data/network/20_gastrointestinal_system.txt/gz
                            -o /path/to/output_directory/
                            -prefix COAD
```
