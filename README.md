# MSIPep:The process of identifying immunogenic peptides
Immunogenic peptide identification pipeline based on proteogenomic strategies for tumor antigen discovery
## Overnew
![SCI-F1-2](https://github.com/user-attachments/assets/e3e69796-1b73-45c7-84e7-eeb90e8fb8e0)

# 1.Introduction

MSIPep is based on immunopeptidomics (MS) data for the identification of tumor antigens. MSIPep consists of four modules. 
Module 1: Preprocess raw RNAseq data, including aligning tumor RNAseq data with corresponding normal tissue RNAseq data or reference genomes to detect mutation sites and generate a theoretical mutant peptide library. HLA typing is also performed. Module 2: Preprocess raw MS data and perform *de novo* peptide sequencing using PepNet. Module 3: Construct a personalized reference protein database. Perform database searches on DDA data using MSFragger and Comet, and identify DIA data using SpectraST and DIA-NN. Module 4: Calculate the binding affinity between peptides and MHC molecules using NetMHCpan, the binding affinity between pMHC molecules and TCR using ProTCR, and assess immunogenicity using DeepImmuno and the IEDB ImmunoGenicity tool, to screen for peptides with high immunogenicity.

# 2. Running environment
MSIPep requires a Linux operation system(centos7) with Python(V3.7), Perl(V5.26) and Java(V1.7) installed.In addition, it is recommended to run the following commands in the terminal to ensure they work correctly.
```
which wget
which pip
which gunzip
```
Before running, make sure the following packages are successfully installed, with the correct versions and mutual compatibility. If you encounter version incompatibility issues between packages, it is recommended to run the scripts in separate Conda environments.

``` 
#PepNet：
Python >= 3.7
Tensorflow >= 2.5.0
Pandas >= 0.20
pyteomics
numba

#Optitype：
1. NumPy
2. Pyomo
3. PyTables
4. Pandas
5. Pysam
6. Matplotlib
7. Future
```
Due to licensing restrictions, **SpectraST** and **DIA-NN** must be obtained separately by the user. **Only when immunopeptidomics data are of the DIA type**, the execution paths of SpectraST and DIA-NN should be specified in `database_search.py`.


# 3. Installation
Before getting started, please run the following code first:
```
chmod -R 777 /path/to/MSIPep/
cd MSIPep
bash start.sh
```

Please place all RNA-seq data for each sample into a separate folder, and place the immunopeptidomics MS data for each sample into a separate folder as well. Make sure all files have sufficient execution permissions.

The population mutation frequency table for different regions is shown below. Please select the corresponding col values based on the region where your sample is located, and use them to specify the col1 and col2 threshold values in the Module 1 run command.
Database name	Crowd	Serial number
| Database name          | Crowd                       | Serial number |
|------------------------|-----------------------------|----------------|
| gnomAD_genome_ALL      | All individuals             | 1              |
| gnomAD_genome_AFR      | African/African American    | 2              |
| gnomAD_genome_AMR      | Latino/Admixed American     | 3              |
| gnomAD_genome_ASJ      | Ashkenazi Jewish            | 4              |
| gnomAD_genome_EAS      | East Asian                  | 5              |
| gnomAD_genome_FIN      | Finnish                     | 6              |
| gnomAD_genome_NFE      | Non-Finnish European        | 7              |
| gnomAD_genome_OTH      | Other                       | 8              |
| ALL.sites.2015_08      | All populations             | 9              |
| EAS.sites.2015_08      | East Asian                  | 10             |
| AFR.sites.2015_08      | African                     | 11             |
| AMR.sites.2015_08      | Admixed American            | 12             |
| EUR.sites.2015_08      | European                    | 13             |
| SAS.sites.2015_08      | South Asian                 | 14             |


# 4.Usage
## 4.1 Module 1 : Preprocessing of RNA-seq data and detection of mutation sites; HLA typing analysis
``` 
python rna_mut_pep.py -1 /patn/to/sample_1.fastq -2 /path/to/sample_2.fastq -t 16 --filter_col1 col1 --filter_col2 col2 --threshold 0.05
#or
python run_mut_pep.py -1 /path/to/sample_1.fastq -2 /path/to/sample_2.fastq -n1 /path/to/normal_1.fastq -n2 /path/to/normal_2.fastq -t 16 --filter_col1 12 --filter_col2 13 --threshold 0.05
#eg.
python rna_mut_pep.py -1 sample_R1.fastq -2 sample_R2.fastq -t 16  --filter_col1 12  --filter_col2 13  --threshold 0.05
``` 
Output the results of RNA-seq data preprocessing and variant calling. If RNA-seq data from the patient’s normal tissue is available, please specify the parameters "-n1" and "-n2" to denote the respective paired-end sequencing reads; the pipeline will then perform mutation detection utilizing both normal and tumor datasets. In the absence of normal tissue data, omit the "-n1" and "-n2" parameters, and the script will default to using the reference genome for variant identification. Additionally, select the appropriate annotation database and designate the corresponding column (col) value from the aforementioned table. The filtering threshold must be calibrated based on the cancer incidence rate of the local population.

## 4.2 Module 2 : Preprocessing of mass spectrometry (MS) data and *de novo* sequencing
``` 
python denovo-pep.py  --input_dir /path/to/raw_or_mgf_files  --output_dir /path/to/output_folder  --output_fasta /path/to/output_folder/final_filtered_peptides.fasta
#eg.
python donovo-pep.py  --input_dir raw_or_mgf_files  --output_dir output_folder  --output_fasta denovo_result.fasta
``` 
If the input data is in RAW format, convert the data format first; if it is in MGF format, perform charge state checking directly. Please ensure that PepNet is functioning properly; if necessary, create a new conda environment.

## 4.3 Module 3 : Database search for peptide identification
```
python database_search.py --fasta_dir1 /path/to/fasta_group1  --fasta_list2 /path/to/extra1.fasta  --mgf /path/to/input.mgf  --output_dir /path/to/output_folder
#eg.
python database_search.py --fasta_dir1 fasta_group1  --fasta_list2 extra.fasta  --mgf input.mgf  --output_dir output_folder
```
The identification results from the database search software are located in the original data folder.

## 4.4 Module 4 : Peptide quantification and filtering
``` 
python immunopep-filter.py /path/to/input_folder /path/to/hla_result_folder
#eg.
python immunopep-filter.py input_folder hla_result_folder
```

The results of high-immunogenicity peptides are stored in the " pep_result " folder.

All required software, reference data files, and mutation annotation files can be downloaded from our Zenodo repository: https://doi.org/10.5281/zenodo.15960309. After downloading, place the files in the specified directory.

# 5.Docker

## 5.1 Pull the Image
You can pull the pre-configured image containing all necessary environments and scripts directly from Docker Hub:
```
docker pull liupeng311/neoantigen-pipeline:msipep
```
## 5.2 Environment Requirements
Docker: Version 20.10.0 or higher is recommended.
Hardware Resources: Due to the intensive nature of variant calling, a minimum of 16 CPU cores and 64GB RAM is suggested.
Data Volume: Always use the -v flag to mount your local data directory to the container.

## 5.3 Quick Start (Usage)
MSIPep consists of several functional modules. You can run specific scripts based on your research needs.

### 5.3.1 RNA Variant Calling & Peptide Generation

Tumor-only mode:
```
docker run --rm -v /data:/data liupeng311/neoantigen-pipeline:msipep \
python rna_mut_pep.py -1 /data/sample_R1.fastq -2 /data/sample_R2.fastq -t 16 --filter_col1 12 --filter_col2 13 --threshold 0.05
```
Tumor-Normal matched mode (Recommended):
```
docker run --rm -v /data:/data liupeng311/neoantigen-pipeline:msipep \
python rna_mut_pep.py -1 /data/tumor_R1.fastq  -2 /data/tumor_R2.fastq -n1 /data/normal_R1.fastq -n2 /data/normal_R2.fastq -t 16 --filter_col1 12 --filter_col2 13  --threshold 0.05
```

### 5.3.2 De novo Peptide Sequencing
```
docker run --rm -v /data:/data liupeng311/neoantigen-pipeline:msipep \
python denovo-pep.py --input_dir /data/raw_or_mgf_files --output_dir /data/output_folder --output_fasta /data/output_folder/final_filtered_peptides.fasta
```

### 5.3.3 Database Search & Validation
```
docker run --rm -v /data:/data liupeng311/neoantigen-pipeline:msipep \
python database_search.py --fasta_dir1 /data/fasta_group1 --fasta_list2 /data/extra1.fasta --mgf /data/input.mgf --output_dir /data/output_folder
```

### 5.3.4 Immunopeptide Filtering
```
docker run --rm -v /data:/data liupeng311/neoantigen-pipeline:msipep \
python immunopep-filter.py data/input_folder /data/hla_result_folder
```
