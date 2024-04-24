# Vibrio-cholerae-PAM


April 23rd 2024

***Full analysis of PAM counts from Vibrio cholerae conjugation experiment 4_17_24.py***

This Python script takes data from eight .fastq files 
Required input files are listed below. These are available from EBI ArrayExpress using the accession code E-MTAB-14044.

The script runs in Python v3 and requires that the following modules are installed:
tkinter
collections
matplotlib.pyplot
numpy

The user is prompted to select eight input files using the tkinter module. The user is also prompted to provide the names and paths of two output files. A third output file is made using a similar path to one of the first two.


------------------------------------------------------------------------------------------------------



List of required input files, with names of example files

1. .fastq for the spacer #4 protospacer library conjugated into a CRISPR-inactive strain, replicate 1

2. .fastq for the spacer #4 protospacer library conjugated into a CRISPR-active strain, replicate 1

3. .fastq for the spacer #21 protospacer library conjugated into a CRISPR-inactive strain, replicate 1

4. .fastq for the spacer #21 protospacer library conjugated into a CRISPR-active strain, replicate 1

5. .fastq for the spacer #4 protospacer library conjugated into a CRISPR-inactive strain, replicate 2

6. .fastq for the spacer #4 protospacer library conjugated into a CRISPR-active strain, replicate 2

7. .fastq for the spacer #21 protospacer library conjugated into a CRISPR-inactive strain, replicate 2

8. .fastq for the spacer #21 protospacer library conjugated into a CRISPR-active strain, replicate 2



------------------------------------------------------------------------------------------------------



List of output files:

1. Tab-delimited .txt file with raw sequence read counts for each of the 65 protospacers from each of the eight samples.

2. Tab-delimited .txt file with normalized conjugation efficiencies for each of the 64 NNN sequences for protospacer #4 and protospacer #21.

3. .png file showing a scatterplot of normalized conjugation efficiencies for each of the 64 NNN sequences for protospacer #4 and protospacer #21.



******************************************************************************************************



------------------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------------

