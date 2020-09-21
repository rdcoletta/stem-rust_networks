# Co-expression Network Analysis of Susceptibility to Wheat Stem Rust

by Rafael Della Coletta and Cory Hirsch (September 2020).

---

The objective of this project is to create gene co-expression networks from wheat and Brachypodium genes in response to wheat stem rust infection.



## Requirements

| Software | Version | Additional libraries / modules |
| -------- | ------- | ------------------------------ |
| R        |         |                                |
| Python   |         |                                |
| Camoco   |         |                                |



## Data

All data, scripts, and output of analyses are located on the folder `/home/hirschc1/della028/projects/stem-rust_networks` from my account at the Minnesota Supercomputing Institute (MSI).

```bash
cd ~/projects/

mkdir -p stem-rust_networks/{analysis,data,scripts}
mkdir -p stem-rust_networks/data/{annotation,expression,go}

# also create a folder to keep MSI output
mkdir -p analysis/msi_dump
```

The lead author of the paper, Eva Henningsen, shared with me all data necessary for building the networks:

`data/annotation`:
* Wheat annotation file (`IWGSC_v1.1_HC_20170706.gff3`)
* Brachypodium annotation file (`BdistachyonBd21_3_537_v1.2.gene_exons.gff3`)

`data/expression`:
* Normalized and raw gene expression matrices for wheat and Brachy (`EH_brachy_counts_normalized.txt`, `EH_wheat_counts_normalized.txt`, `wheat_counts_raw.txt`, `brachy_counts_raw.txt`)

`data/go`:
* GO terms for wheat and Brachy (`wheat_full_gos_long.txt`, `brachy_full_gos_long.txt`)



### Normalize raw expression by FPKM

Camoco requires that expression data is normalized by FPKM, but the normalized files Eva sent to me was generated by DEseq2 package from R, which is not FPKM. Thus, I wrote the script `scripts/counts2fpkm.R` to normalize the raw counts from HTseq into FPKM according to guidelines described by the [National Cancer Institute Genomic Data Commons](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#mrna-expression-ht-seq-normalization). It uses the following formula to calculate FPKM:

```
                     (reads mapped to the gene) * 10^9
 FPKM = ------------------------------------------------------------
         (reads mapped to all protein-coding genes) * (gene length)
```

I ran `scripts/counts2fpkm.R` for both wheat and brachy datasets.

```bash
# first need to remove header from annotation
grep -v "^#" data/annotation/IWGSC_v1.1_HC_20170706.gff3 > data/annotation/IWGSC_v1.1_HC_20170706.no-header.gff3
grep -v "^#" data/annotation/BdistachyonBd21_3_537_v1.2.gene_exons.gff3 > data/annotation/BdistachyonBd21_3_537_v1.2.gene_exons.no-header.gff3

# then remove quotes around gene and sample names
sed -i 's/"//g' data/expression/wheat_counts_raw.txt
sed -i 's/"//g' data/expression/brachy_counts_raw.txt
# remove htseq-specific lines (i.e. start with _)
sed -i '/^_/d' data/expression/wheat_counts_raw.txt
sed -i '/^_/d' data/expression/brachy_counts_raw.txt

# for how to use script
Rscript scripts/counts2fpkm.R --help

# wheat dataset
Rscript scripts/counts2fpkm.R data/expression/wheat_counts_raw.txt \
                              data/annotation/IWGSC_v1.1_HC_20170706.no-header.gff3 \
                              data/expression/wheat_counts_fpkm.txt \
                              --cores=10

# brachy dataset
Rscript scripts/counts2fpkm.R data/expression/brachy_counts_raw.txt \
                              data/annotation/BdistachyonBd21_3_537_v1.2.gene_exons.no-header.gff3 \
                              data/expression/brachy_counts_fpkm.txt \
                              --cores=10
```

As a quick QC, I wrote `scripts/pca_expr_data.R` to perform PCA on the expression datasets (raw counts and fpkm) for wheat and Brachy.

```bash
mkdir -p analysis/qc/expression

for geno in wheat brachy; do
  for expr in raw fpkm; do
    Rscript scripts/pca_expr_data.R data/expression/${geno}_counts_${expr}.txt \
                                    analysis/qc/expression/pca_${geno}_${expr}.png
  done
done
```

> Grouping of samples according to PCA seems to agree reasonably well.



### Filter data by CV and samples for each network

Coexpression networks require variation of gene expression across samples. Thus, filtering genes by their coefficient of variation (CV) eliminates such uninformative genes, reduces the initial number of genes to build networks and speeds up clustering. I also have to filter the whole dataset to have the samples I want for each network. There will be **one network for each genotype** (wheat susceptible, wheat resistant, brachy) using all reps for 2, 4 and 6 days after inoculating the plants with stem rust, and all reps for the 2dpi samples that were mock inoculated.

To have an idea of what CV threshold to use, I plotted the distribution of genes CV among samples using `scripts/prepare_data_for_camoco.R`.

```bash
# plot CV distribution among all samples
Rscript scripts/prepare_data_for_camoco.R data/expression/brachy_counts_fpkm.txt \
                                          analysis/qc/expression/cv_distribution_brachy-fpkm.png \
                                          --plot-cv
Rscript scripts/prepare_data_for_camoco.R data/expression/wheat_counts_fpkm.txt \
                                          analysis/qc/expression/cv_distribution_wheat-fpkm.png \
                                          --plot-cv

# add sample names to variables
brachy_samples="Bd21_D2_mock_R1,Bd21_D2_mock_R2,Bd21_D2_mock_R3,Bd21_D2_treated_R1,Bd21_D2_treated_R2,Bd21_D2_treated_R3,Bd21_D4_treated_R1,Bd21_D4_treated_R2,Bd21_D4_treated_R3,Bd21_D6_treated_R1,Bd21_D6_treated_R2,Bd21_D6_treated_R3"
wheatR_samples="Sr9b_D2_mock_R1,Sr9b_D2_mock_R2,Sr9b_D2_mock_R3,Sr9b_D2_treated_R1,Sr9b_D2_treated_R2,Sr9b_D2_treated_R3,Sr9b_D4_treated_R1,Sr9b_D4_treated_R2,Sr9b_D4_treated_R3,Sr9b_D6_treated_R1,Sr9b_D6_treated_R2,Sr9b_D6_treated_R3"
wheatS_samples="W2691_D2_mock_R1,W2691_D2_mock_R2,W2691_D2_mock_R3,W2691_D2_treated_R1,W2691_D2_treated_R2,W2691_D2_treated_R3,W2691_D4_treated_R1,W2691_D4_treated_R2,W2691_D4_treated_R3,W2691_D6_treated_R1,W2691_D6_treated_R2,W2691_D6_treated_R3"

# plot CV distribution by genotype
Rscript scripts/prepare_data_for_camoco.R data/expression/brachy_counts_fpkm.txt \
                                          analysis/qc/expression/cv_distribution_brachy-fpkm.inf_2-4-6.mock_2.png \
                                          --plot-cv \
                                          --keep-samples=$brachy_samples
Rscript scripts/prepare_data_for_camoco.R data/expression/wheat_counts_fpkm.txt \
                                          analysis/qc/expression/cv_distribution_Sr9b-fpkm.inf_2-4-6.mock_2.png \
                                          --plot-cv \
                                          --keep-samples=$wheatR_samples
Rscript scripts/prepare_data_for_camoco.R data/expression/wheat_counts_fpkm.txt \
                                          analysis/qc/expression/cv_distribution_W2691-fpkm.inf_2-4-6.mock_2.png \
                                          --plot-cv \
                                          --keep-samples=$wheatS_samples
```

Below is the number of genes that would be eliminated by different CV thresholds for all samples. Based on this table and the distribution of CV for each genotype, I will only keep genes with CV >= 0.1 for both brachy and wheat.

|               | Brachy               | Wheat                |
| ------------- | -------------------- | -------------------- |
| Not expressed | 8398 genes (21.5%)   | 28795 genes (26.69%) |
| CV < 0.1      | 8518 genes (21.8%)   | 28797 genes (26.69%) |
| CV < 0.2      | 12678 genes (32.45%) | 30307 genes (28.09%) |
| CV < 0.3      | 18525 genes (47.42%) | 38032 genes (35.25%) |
| CV < 0.4      | 22511 genes (57.62%) | 47331 genes (43.87%) |
| CV < 0.5      | 25351 genes (64.89%) | 55537 genes (51.48%) |

An additional formatting of the expression matrix is required for building networks with Camoco. `scripts/prepare_data_for_camoco.R` generates a `.csv` file without the first comma from the header, in order to the input dataset look like this:

| Sample1 | Sample2 | Sample3 | Sample4 |    |
|---------|---------|---------|---------|----|
| gene1   | 10      | 20      | 30      | 40 |
| gene2   | 10      | 20      | 30      | 40 |
| gene3   | 10      | 20      | 30      | 40 |

```bash
# brachy dataset
Rscript scripts/prepare_data_for_camoco.R data/expression/brachy_counts_fpkm.txt \
                                          data/expression/expr_data_fpkm_brachy.cv_0.1.inf_2-4-6.mock_2.csv \
                                          --filter-cv=0.1 \
                                          --keep-samples=$brachy_samples

# wheat resistant dataset
Rscript scripts/prepare_data_for_camoco.R data/expression/wheat_counts_fpkm.txt \
                                          data/expression/expr_data_fpkm_Sr9b.cv_0.1.inf_2-4-6.mock_2.csv \
                                          --filter-cv=0.1 \
                                          --keep-samples=$wheatR_samples

# wheat susceptible dataset
Rscript scripts/prepare_data_for_camoco.R data/expression/wheat_counts_fpkm.txt \
                                          data/expression/expr_data_fpkm_W2691.cv_0.1.inf_2-4-6.mock_2.csv \
                                          --filter-cv=0.1 \
                                          --keep-samples=$wheatS_samples
```



### Adjust format of Brachy gene IDs

After some preliminary tests, I found out that the gene IDs from Brachy genome annotation is not compatible with Camoco (it couldn't match the gene IDs of the reference annotation with the IDs of the expression data or Gene Ontology), very likely due to the presence of non-alphabetic characters such as `-` and `.`. I realized that eliminating these characters from reference annotation file, Gene Ontology file, and gene expression data solved the problem. For example, gene ID `BdiBd21-3.1G0000100.v1.2` became `BdiBd1G0000100`.

```bash
# reference genome annotation
grep -P "^#" data/annotation/BdistachyonBd21_3_537_v1.2.gene_exons.gff3 > data/annotation/BdistachyonBd21_3_537_v1.2.gene_exons.corrected-IDs.gff3
awk 'BEGIN{FS=OFS="\t"} {gsub(/\.v1\.2/, "", $9)} 1' data/annotation/BdistachyonBd21_3_537_v1.2.gene_exons.no-header.gff3 | awk 'BEGIN{FS=OFS="\t"} {gsub(/BdiBd21-3\./, "BdiBd", $9)} 1' >> data/annotation/BdistachyonBd21_3_537_v1.2.gene_exons.corrected-IDs.gff3

# awk commands:
#   BEGIN{} this block of code will be executed before processing any input line
#   FS=OFS="\t" set input and output field separator as tab
#   gsub(/\.v1\.2/, "", $9 for each input line, replace all the '.v1.2' in 9th field with nothing
#   1 is an awk idiom to print contents of $0 (which contains the input record)

# go file
sed 1d data/go/brachy_full_gos_long.txt | awk 'BEGIN{FS=OFS="\t"} {gsub(/BdiBd21-3\./, "BdiBd", $1)} 1' > data/go/brachy_full_gos_long.corrected-IDs.txt

# expression data
head -n 1 data/expression/expr_data_fpkm_brachy.cv_0.1.inf_2-4-6.mock_2.csv > data/expression/expr_data_fpkm_brachy.cv_0.1.inf_2-4-6.mock_2.corrected-IDs.csv
sed 1d data/expression/expr_data_fpkm_brachy.cv_0.1.inf_2-4-6.mock_2.csv | awk 'BEGIN{FS=OFS="\t"} {gsub(/\.v1\.2/, "", $1)} 1' | awk 'BEGIN{FS=OFS="\t"} {gsub(/BdiBd21-3\./, "BdiBd", $1)} 1' >> data/expression/expr_data_fpkm_brachy.cv_0.1.inf_2-4-6.mock_2.corrected-IDs.csv
```


## Install Camoco

Camoco was already installed in a **virtual environment** using the software [miniconda](https://conda.io/miniconda.html).

```bash
cd ~/software/

# get miniconda bash installer
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
# install miniconda
bash Miniconda3-latest-Linux-x86_64.sh -p ~/software/miniconda
# check if conda needs uptdates
conda update -n base conda
# create virtual environment for camoco
conda create -n camoco python=3.6
# check installed environments
conda info --envs

# update pip, wheel and setuptools
pip install --upgrade pip
pip install --upgrade wheel
pip install --upgrade setuptools

# activate environment and install Camoco on that environment
source activate camoco
# notice that (camoco) now appears before "della028@ln0005 [~/software] %"
pip install numpy
pip install camoco
# test if package was installed correctly
camoco -h

# deactivate environment when done working
source deactivate
```



### Load reference annotation and gene ontology

Before building any network, I need to load Brachy and wheat reference genome annotations into Camoco's database. Once a reference genome is loaded, you don't need to load it again: it will stay in Camoco's database to be used as many time as you want.


```bash
# activate camoco virtual environment
source activate camoco

# create a reference genome dataset for camoco -- see 'camoco build-refgen -h' for help
camoco build-refgen data/annotation/BdistachyonBd21_3_537_v1.2.gene_exons.corrected-IDs.gff3 BrachyRef Bd21_3_537_v1.2 phytozomev13 Brachypodium_distachyon
camoco build-refgen data/annotation/IWGSC_v1.1_HC_20170706.gff3 WheatRef TraesChineseSpring_v1.1_201706 IWGSC Triticum_aestivum
# check that loading was successfull
camoco ls

# deactivate virtual environment
source deactivate
```

Camoco performs Gene Ontology (GO) enrichment analysis to quality control the networks created. In order to do that, I need to create a GO object in Camoco, which requires a `go.obo` file containing all core ontology terms (<http://geneontology.org/docs/download-ontology/>) and a two-column species-specific file relating genes (first column) and their respective GO terms (second column).

The `go.obo` file is a comprehensive list of GO terms, thus it contains very complex and specific terms that may be hard to analyze or undertand what's really going on. Thus, one alternative is to use more generic terms (`goslim_generic.obo`) or use terms specific to plants (`goslim_plant.obo`; which is derived from Arabidopsis).

```bash
# activate camoco virtual environment
source activate camoco

# get obo files
wget -P data/go/ http://purl.obolibrary.org/obo/go.obo
wget -P data/go/ http://current.geneontology.org/ontology/subsets/goslim_generic.obo
wget -P data/go/ http://current.geneontology.org/ontology/subsets/goslim_plant.obo

# remove header from wheat go files
sed -i 1d data/go/wheat_full_gos_long.txt

# create a go annotation reference for camoco  -- see 'camoco build-go -h' for help
camoco build-go data/go/brachy_full_gos_long.corrected-IDs.txt data/go/go.obo BrachyGO brachy_go_annotation BrachyRef
camoco build-go data/go/wheat_full_gos_long.txt data/go/go.obo WheatGO wheat_go_annotation WheatRef

# also create generic go database
camoco build-go data/go/brachy_full_gos_long.corrected-IDs.txt data/go/goslim_generic.obo BrachyGOslim brachy_go_generic BrachyRef
camoco build-go data/go/wheat_full_gos_long.txt data/go/goslim_generic.obo WheatGOslim wheat_go_generic WheatRef

camoco build-go data/go/brachy_full_gos_long.corrected-IDs.txt data/go/goslim_plant.obo BrachyGOplant brachy_go_plant BrachyRef
camoco build-go data/go/wheat_full_gos_long.txt data/go/goslim_plant.obo WheatGOplant wheat_go_plant WheatRef


# check that loading was successfull
camoco ls

# deactivate virtual environment
source deactivate
```

> Done: Ontology:BrachyGO - desc: brachy_go_annotation - contains 3589 terms for Reference Genome: Brachypodium_distachyon - phytozomev13 - BrachyRef
> Done: Ontology:WheatGO - desc: wheat_go_annotation - contains 12572 terms for Reference Genome: Triticum_aestivum - IWGSC - WheatRef





## Build networks

Originally, Meesh had generated three networks, one for each genotype, and included infected samples 2,4 and 6 dpi and mock samples 2dpi only. To build the networks with `scripts/build_network.sh`, I need to specify five parameters: the Camoco filters to be applied to expression dataset (`OPT`), the expression data file (`IN`), a name for the network (`NAME`), a description of the network (`DESC`) and the reference genome to be used (`REF`). An additional parameter (`HEALTH`) is a folder to save summary statistics of the network.

```bash
mkdir -p analysis/qc/health

# set parameters build networks
options="--max-gene-missing-data 0.4 --max-accession-missing-data 0.4 --min-single-sample-expr 0.5 --min-expr 0.001"

# brachy
qsub -v OPT="$options",IN=data/expression/expr_data_fpkm_brachy.cv_0.1.inf_2-4-6.mock_2.corrected-IDs.csv,NAME=SR_brachy_1,DESC=brachy_infected-2-4-6dpi_mock_2dpi_mgmd0-4_mamd0-4_msse0-5_min0-001,REF=BrachyRef,HEALTH=analysis/qc/health scripts/build_network.sh
sleep 60
# wheat R
qsub -v OPT="$options",IN=data/expression/expr_data_fpkm_Sr9b.cv_0.1.inf_2-4-6.mock_2.csv,NAME=SR_Sr9b_1,DESC=Sr9b_infected-2-4-6dpi_mock_2dpi_mgmd0-4_mamd0-4_msse0-5_min0-001,REF=WheatRef,HEALTH=analysis/qc/health scripts/build_network.sh
sleep 60
# wheat S
qsub -v OPT="$options",IN=data/expression/expr_data_fpkm_W2691.cv_0.1.inf_2-4-6.mock_2.csv,NAME=SR_W2691_1,DESC=W2691_infected-2-4-6dpi_mock_2dpi_mgmd0-4_mamd0-4_msse0-5_min0-001,REF=WheatRef,HEALTH=analysis/qc/health scripts/build_network.sh
```

> NOTE 1: Use `OPT="$options"` instead of `OPT=$options` because of the whitespaces inside `$options`.
> NOTE 2: Starting multiple jobs at the same time will cause an error due to accessing Camoco database at the same time, that's why I added `sleep 60` between `qsub` commands.

Camoco can also perform GO enrichment analysis of the entire network as QC. To do so, I wrote `scripts/network_health.sh`.

```bash
# # summary with GO enrichment analysis -- take much longer!
# qsub -v NAME=SR_brachy_1,REF=BrachyRef,GO=BrachyGO,OUT=analysis/qc/health/SR_brachy_1 scripts/network_health.sh
# qsub -v NAME=SR_Sr9b_1,REF=WheatRef,GO=WheatGO,OUT=analysis/qc/health/SR_Sr9b_1 scripts/network_health.sh
# qsub -v NAME=SR_W2691_1,REF=WheatRef,GO=WheatGO,OUT=analysis/qc/health/SR_W2691_1 scripts/network_health.sh
```

After networks are built, I retrieved cluster information (i.e. which genes are in which clusters) from Camoco database with `scripts/retrieve_network_info.py`. The cluster number is sorted by cluster size, meaning that cluster 0 is the largest cluster, cluster 1 is the second largest and so on. In addition, this script also retrieves the transformed expression data used as input for building networks (Camoco performs an inverse hyperbolic sine transformation due to the dynamic range of the RNAseq data).

```bash
mkdir -p analysis/clusters

source activate camoco

for network in SR_brachy_1 SR_Sr9b_1 SR_W2691_1; do
  python scripts/retrieve_network_info.py $network analysis/clusters
done

source deactivate
```
