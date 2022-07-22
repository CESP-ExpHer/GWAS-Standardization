# GWAS-Standardization
When you work on **GWAS Summary Statistics data** in **more than one trait**, it is important to put alleles in both data in a same position. We called this step as **Standardization**. <br>
Here, we provided two Python script in order to clean, order, and standardize your data:
## First Trait
The code reads the GWAS Summary Statistics data for **BCAC_2020_onco_all** *"icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt"* and writes it in a proper format that could be used for our future analyses. 
These are headers of the input file that we use in the script:
<br>
| SNP_ID | chr | position | Effect_Allele | non-Effect_Allele | eaf | info | beta | se | P_value |
| -- | -- | -- | -- | -- | -- | -- | -- | -- | -- |

Here are the steps in the file:
- Removing duplicated SNPs

### libraries used in the code
```python
from collections import defaultdict
from scipy import stats
import datetime
```
## Second Trait
