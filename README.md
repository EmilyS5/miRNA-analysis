# miRNA-analysis
Program to analyse miRNA-gene interactions and identify novel diagnostic biomarkers

- Bioinformatics pipeline integrating three databases (VarElect, miRDB, and miRDIP) to analyse miRNA-gene interactions and identify diagnostic biomarkers and their regulatory networks.
- This offers a scalable approach to assess miRNA-gene interactions for viral infections
- The use of this program was validated through data from experiments into novel diagnostic biomarkers for HRV16.

To use:
- VarElect GeneCards API required to allow for mapping of the miRNA gene targets to the chosen phenotype (RV in this program)
- miRDB is queried using locally stored data, which can be downloaded from https://mirdb.org/download.html (miRDB v6.0, June 2019, MirTargetV4, miRbase 22)
- Program coded in python (3.11.9)
- Update miRNAs, desired number of gene targets, phenotype, and file names as stated in the code comments.

Program run through python powershell as administrator to allow for the CSV to be saved.

Expected output:
- CSV file containing list of top 10 gene targets per miRNA per database (20 total per miRNA)
- CSV file containing list of gene targets found to be relevant to the phenotype 
