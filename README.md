# GWAS-Standardization
Created by: Yazdan Asgari<br>
Creation date: 22 Jul 2022<br>
Update: Jul 2022<br>
https://cesp.inserm.fr/en/equipe/exposome-and-heredity
<br>
<br>

## Contents
- [First Trait](#first-trait)
  * [INPUT](#input)
  * [SCRIPT](#script)
  * [OUTPUTS](#outputs)
  * [libraries used in the code](#libraries-used-in-the-code)
- [Second Trait](#second-trait)
  * [INPUTS](#inputs)
  * [SCRIPT](#script-1)
  * [OUTPUTS](#outputs-1)
  * [libraries used in the code](#libraries-used-in-the-code-1)

When you work on **GWAS Summary Statistics data** in **more than one trait**, it is important to put alleles in both data in a same position. We called this step as **Standardization**. <br><br>

Here, we provided two Python script in order to clean, order, and standardize your data:
## First Trait
**IMPORTANT NOTE:** For the first trait, we assume the data for Breast cancer consortium ([BCAC version 2020](https://bcac.ccge.medschl.cam.ac.uk/bcacdata/oncoarray/oncoarray-and-combined-summary-result/)).<br><br>
**NOTE:** The code reads the GWAS Summary Statistics data for **BCAC_2020_onco_all** *"icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt"* and writes it in a proper format that could be used for our future analyses. <br><br>

### INPUT
These are headers of the input file that we use in the script:

| SNP_ID | chr | position | Effect_Allele | non-Effect_Allele | eaf | info | beta | se | P_value |
| -- | -- | -- | -- | -- | -- | -- | -- | -- | -- |

eaf = Effect Allele Frequency
<br>
<br>

### SCRIPT
Here are the steps in the file:
- reading the file
- if (eaf or info or beta or se or P_value) = NULL or NA, it removes the SNP (means it skips the line)
- if both beta and se = 0, it removes the SNP (means it skips the line)
- if one of the alleles has a length not equal to 1, and the SNP includes alleles that are not composed of ATGC, it will be removed 
- if the SNP is an ambiguous SNPs, it removes the SNP (means it skips the line) and prints the SNP in its output file
- if A1 and A2 have more than one nucleotide and have the same length, it checks if they are not ambiguous. If so, it removes the SNP (means it skips the line)
- if the SNP is in the duplicated set, it removes the SNP (means it skips the line) and prints the SNP in its output file
- it calculates the MAF (Minor Allele Frequency) based on EAF 
- it counts the SNPs with a value smaller than 0.9 for info

### OUTPUTS
Output file: writing the line in the output file which includes the following headers names:

| ID | CHR | POS | A1 | A2 | beta | se | p | info | ngt | eaf | MAF |
| -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- |

A1 = Effect Allele
<br>
eaf = Effect Allele Frequency
<br>
A2 = non-Effect Allele
<br>
ngt = used for LD score tool
<br>
<br>
The script keeps track of the different categories of SNPs in an output file:
- SNPs which kept, 
- SNPs with a low value for info (or r2),
- SNPs which are ambiguous, 
- SNPs with NULL values for beta/se/info/pval,
- SNPs with zero beta&se,
- SNPs which are duplicated,# SNPs include an allele that are not composed of ATGC,
- SNPs with an allele that are not composed of ATGC,
- SNPs with a length > 1 and with alleles that are not composed of ATGC

### libraries used in the code
```python
from collections import defaultdict
from scipy import stats
import datetime
```

**NOTE:** The code (*A1_code_reformatting_file_bcac_2020_all.py*) is available here ([Download](https://github.com/CESP-ExpHer/GCPBayes-Pipeline/blob/main/0_Codes))

## Second Trait

**IMPORTANT NOTE:** For the second trait, we assume the data for Ovarian cancer consortium ([OCAC version 2020](https://ocac.ccge.medschl.cam.ac.uk/)).
<br><br>
**NOTE:** The code reads the GWAS Summary Statistics data for OCAC *"extraction_OCAC.txt"* and writes it in a proper format that could be used for our future analyses. It also checks the *"Effect"* and *"non-Effect"* Alleles to be in accordance with the **reference file** (here the file created from above section (**First Trait**))
<br><br>

### INPUTS
1. The file created from above section (**For the First Trait**)
2. The file that needs Standardization (**For the Second Trait**)
<br>
These are headers of the input file that needs Standardization:

| chr | position | A1 | A2 | EAF | nEAF | beta | se | p-value | info | N |
| -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- |

A1 = Effect Allele
<br>
A2 = non-Effect Allele
<br>
EAF = Effect Allele Frequency
<br>
nEAF = non-Effect Allele Frequency
<br>
N = Number of samples
<br>
<br>

### SCRIPT
Here are the steps in the file:
- reading the file
- if (EAF or nEAF or beta or se or p) = NULL or NA, it removes the SNP (means it skips the line)
- if both beta and se = 0, we skip the line
- if the SNP is in the duplicated set, it removes the SNP (means it skips the line) and prints the SNP in its output file
- if A1 and A2 have more than one nucleotide and have the same length, it checks if they are not ambiguous. If so, it removes the SNP (means it skips the line)
- Standardization (Alleles Checking) Steps
  - checking the reference file to see if it exists (means if CHR:POS exists in the reference file)
  - if one of the alleles has a length not equal to 1, and the SNP includes alleles that are not composed of ATGC, it will be removed
  - if the SNP is an ambiguous SNPs, it removes the SNP (means it skips the line) and prints the SNP in its output file
  - if the alleles of the file and the reference alleles are the same, storing the A1 and A2 from the reference file, and writing the line in the output file
  - if A1 and A2 are inverse from the reference, storing the A1 and A2 from the reference file in the right order. Also, it changes the sign of the beta and writing the line in the output file
  - if the alleles of the file are the complements of the reference alleles, storing the A1 and A2 from the reference file, and writing the line in the output file
  - if the alleles of the file are the complements of the reference alleles, BUT they are switched, storing the A1 and A2 from the reference file. Also, it changes the sign of the beta and writing the line in the output file
  - If the SNP does not fit to any of the above categories, it removes the SNP (means it skips the line)
- if the SNP does not exist in the reference, it removes the SNP (means it skips the line)

### OUTPUTS
Output file: writing the line in the output file which includes the following headers names:

| snp | chr | bp_hg19 | Effect_A | nonEffect_A | beta | se | pval | info | ngt | EAF | nEAF | MAF |
| -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- |

EAF = Effect Allele Frequency
<br>
nEAF = non-Effect Allele Frequency
<br>
ngt = used for LD score tool
<br>
<br>

The script keeps track of the different categories of SNPs in an output file:
- SNPs which kept, 
- SNPs which are ambiguous, 
- SNPs with NULL values for freq/beta/se/pval,
- SNPs with zero beta&se,
- SNPs which are duplicated,
- SNPs include an allele that are not composed of ATGC,
- SNPs with allele size more than one, 
- SNPs which removed due to mismatches of both alleles between the file and reference files, 
- SNPs which removed due to not be available in the reference file

### libraries used in the code
```python
from collections import defaultdict
import datetime
```

**NOTE:** The code (*A2_code_reformatting_file_ocac_bcac_2020_all.py*) is available here ([Download](https://github.com/CESP-ExpHer/GCPBayes-Pipeline/blob/main/0_Codes))

## Checking The Created Files
In order to make sure if the scripts worked well and the files are correct (in terms of alleles and values), you could use **"CheckSumStats"**. There is a very detailed tutorial for how to use it in their [GitHub page](https://github.com/MRCIEU/CheckSumStats).


