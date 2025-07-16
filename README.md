# MSIPep:The process of identifying immunogenic peptides
Immunogenic peptide identification pipeline based on proteogenomic strategies for tumor antigen discovery
## Overnew
![SCI-F1-2](https://github.com/user-attachments/assets/e3e69796-1b73-45c7-84e7-eeb90e8fb8e0)

MSIPep, based on immunopeptide MS data. By analyzing the RNAseq data and immunopeptide MS data of patients, tumor neoantigens from the coding and non-coding regions were identified, including peptide segments that might come from circRNA and lncRNA. MSIPep consists of four modules. Module 1: Preprocessing of raw RNAseq data, which involves sequence alignment of tumor RNAseq data with the corresponding normal tissue RNAseq data or reference genomes to detect mutation sites and generate a theoretical mutation peptide library. And perform HLA typing calculation; Module 2: The original MS data is preprocessed and de novo sequencing is performed using PepNet. The peptide segments are preliminarily filtered through the two parameters of Score and PPM Difference. Module 3: Build a personalized reference protein database, conduct database search on DDA data using MAFragger and Comet, and identify DIA data using SpectraST and DIA-NN. Module 4: The binding affinity between peptides and MHC molecules was calculated using NetMHCpan, the binding affinity between pMHC molecules and TCR was calculated using ProTCR, and the Immunogenicity was calculated using DeepImmuno and IEDB tools ImmunoGenicity. Peptides with high immunogenicity were screened out. Compared with other processes, MSIPep has a lower threshold for input data and can provide more comprehensive identification results.

2. Running environment
MSIPep requires a Linux operation system(centos7) with Python(V3.7), Perl(V5.26) and Java(V1.7) installed.In addition, it is recommended to run the following commands in the terminal to ensure they work correctly.
''' 
which wget
which pip
which gunzip
 '''

Before getting started, please run the following code first:
```
cd MSIPep
bash start.sh
```



## Table of Contents
1.Features
2.Install
3.Usage
4.Pipeline Modules
5.Output Structure
6.License
7.Contact
##Features
MSIPep provides

##Install
NetMHCpan 
DIA-NN

##Usage
jiang shu ju an yi xia ge shi fang zai data wen jain jia zhong 

##Pipeline Modules
RNA
MS
Peptide

##Output Structure
shuchushili 

##License

##Contact
