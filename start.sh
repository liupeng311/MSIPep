#!/bin/bash

# Switch to the "software" directory
cd software || exit 1

# Install common software packages (note: the correct spelling is 'install')
conda install -y bwa
conda install -y samtools
conda install -y gatk
conda install -y picard
conda install -y razers3

# Download necessary files
wget -O TPP_7.3.0-src.tgz https://sourceforge.net/projects/sashimi/files/Trans-Proteomic%20Pipeline%20%28TPP%29/TPP%20v7.3%20%28Trade%20Wind%29%20rev%200/TPP_7.3.0-src.tgz/download

wget -O annovar.zip https://zenodo.org/records/15960309/files/annovar.zip?download=1
wget -O comet.zip https://zenodo.org/records/15960309/files/comet.zip?download=1
wget -O DeepImmuno.zip https://zenodo.org/records/15960309/files/DeepImmuno.zip?download=1
wget -O DIA-NN-2.2.0.zip https://zenodo.org/records/15960309/files/DIA-NN-2.2.0.zip?download=1
wget -O humandb.zip https://zenodo.org/records/15960309/files/humandb.zip?download=1
wget -O immunogenicity.zip https://zenodo.org/records/15960309/files/immunogenicity.zip?download=1
wget -O msconvert.zip https://zenodo.org/records/15960309/files/msconvert.zip?download=1
wget -O MSFragger-3.8.zip https://zenodo.org/records/15960309/files/MSFragger-3.8.zip?download=1
wget -O ncbi-blast.zip https://zenodo.org/records/15960309/files/ncbi-blast.zip?download=1
wget -O netMHCpan-4.1.zip https://zenodo.org/records/15960309/files/netMHCpan-4.1.zip?download=1
wget -O OptiType.zip https://zenodo.org/records/15960309/files/OptiType.zip?download=1
wget -O PepNet.zip https://zenodo.org/records/15960309/files/PepNet.zip?download=1
wget -O reference.zip https://zenodo.org/records/15960309/files/reference.zip?download=1
wget -O Trimmomatic.zip https://zenodo.org/records/15960309/files/Trimmomatic.zip?download=1

# Unzip all zip files
unzip -o '*.zip'

# Move the "humandb" and "reference" folders up one level
mv humandb ../
mv reference ../

# Unzip nested zip files inside "reference" directory
cd ../reference || exit 1
unzip -o '*.zip'

# Unzip nested zip files inside "humandb" directory
cd ../humandb || exit 1
unzip -o '*.zip'

#Build reference genome index with BWA
bwa index ./reference/hg38.fa
