# GWAS-Standardization
When you work on **GWAS Summary Statistics data** in **more than one trait**, it is important to put alleles in both data in a same position. We called this step as **Standardization**. <br>
Here, we provided two Python script in order to clean, order, and standardize your data:
## First Trait
**NOTE:** The code reads the GWAS Summary Statistics data for **BCAC_2020_onco_all** *"icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt"* and writes it in a proper format that could be used for our future analyses. 
These are headers of the input file that we use in the script:
<br>
| SNP_ID | chr | position | Effect_Allele | non-Effect_Allele | eaf | info | beta | se | P_value |
| -- | -- | -- | -- | -- | -- | -- | -- | -- | -- |

Here are the steps in the file:
- if (eaf or info or beta or se or P_value) = NULL or NA, it removes the SNP (means it skips the line)
- if both beta and se = 0, it removes the SNP (means it skips the line)
- if one of the alleles has a length not equal to 1, and the SNP includes alleles that are not composed of ATGC, it will be removed 
- if the SNP is an ambiguous SNPs, it removes the SNP (means it skips the line) and prints the SNP in its output file
- if A1 and A2 have more than one nucleotide and have the same length, it checks if they are not ambiguous. If so, it removes the SNP (means it skips the line)
- if the SNP is in the duplicated set, it removes the SNP (means it skips the line) and prints the SNP in its output file
- it calculates the MAF (Minor Allele Frequency) based on EAF 
- it counts the SNPs with a value smaller than 0.9 for info
- writing the line in the output file which includes the following headers names

| ID | CHR | POS | A1(Effect_Allele) | A2(non-Effect_Allele) | beta | se | p | info | 0 | eaf | MAF |
| -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- | -- |

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
