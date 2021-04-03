# GMMchi

GMMchi is the python package for Gaussian Mixture Modeling using the chi-square protocol. GMMchi enables the efficient subcategorization gene expression data in large datasets. The method is based on identifying mixtures of normal and non-normal distributed tails. Although GMMchi is developed to identify unique patterns in gene expression data, the general use case of pattern identification and subcategorization creates a much wider application for GMMchi in any datasets that exhibit mixtures of normal or non-normal data.

When applying GMMchi, a-priori knowledge of distinct subpopulations due to underlying mechanisms (i.e. mutation, methylation..etc) is ideal for interpreting the result of the analysis. The pre-print/soon-to-be published use-case was meant for gene expression data analysis where bimodal distributions are related in relation to mutated v.s. normal populations. In other words, the expectation is that the mutated and wildtype subpopulation expresses levels of expression that cluster in distinct distributions. 

GMMchi offers a systematic approach for identifying and characterizing different patterns of normal and non-normal mixtures. Moreover, the advantage of transforming continuous data into categorized data enables researchers/users the ability to study and correlate genes to phenotypes and to explore data using a pattern-based (vs the traditional trend-base) analysis. This package assumes users with zero python knowledge thus starts with very simple instructions. We will explore several examples of GMMchi on gene expression analysis.

## Getting Started 

Download Anaconda at https://www.anaconda.com/distribution/

After downloading, open terminal (Mac) or cmd (Windows). Open Jupyter Notebook either by entering 'jupyter notebook' in the terminal/cmd or using the anaconda application downloaded. This will open an IDE (Integrated development environment) using your default browser. Jupyter notebook is essentially representing your computer files in an IDE and allows you to directly interact with blocks of python code, making it a much more pleasant experience vs using the cli (command line interface). 

Create a new python script by creating a new folder by clicking new --> python 3 notebook located on the top right corner. You can also use an existing folder where you keep your python scripts. 

### Download Package

Download the GMMchi package by:
```
pip install git+https://github.com/jeffliu6068/GMMchi.git
```
or 
```
pip install GMMchi
```

### Import

Once installed, import the package by: 

```
import GMMchi
```
## Intuition: How GMMchi Works in Gene Expression Analysis

The goal of GMMchi in gene expression analysis is the categorization of continuous data into 1s, 2s and occasionally 3s. 1s and 2/3s represent low or non-expressing vs high-expressing samples of any given gene, respectively. There are many ways GMMchi-categorized data can be studied in downstream analysis. The steps outlined below is an example of a standard method of analyzing a large dataset: 

### Postprocessing

1) Determine the background threshold of your input sample
2) Filter and remove genes that are not expressed by any of the samples

### Analysis

3) If doing analysis on a single gene, categorize a gene expression by applying GMMchi on your gene of interest in the postprocessed data
4) If doing a full-scale screen or analysis, categorize all gene expresssions by running a for loop to apply GMMchi on each gene to recreate a categorized matrix with each sample cateogrized as 1, 2, or 3
5) Run 2x2 table analysis on the categorized matrix 

# Available Tools in the GMMchisquare Package

## Calculating Background Threshold

```
means, std, filt = GMMchi.GMM_modelingt('TCGA Colorectal Cancer', input_data_cancer, log2transform=True, 
                        verbosity = True, Single_tail_validation=False, calc_back = True)
```
### Input

'TCGA Colorectal Cancer': When calc_back = True, input string will be automatically used as the title for your output graphs

input_data_cancer: Dataframe with genes (row) x samples (columns)
             
log2transform: perform log2-transformation on the data

verbosity: print each stage of GMMchi

calc_back: Boolean to indicate whether to calculate background threshold

Single_tail_validation: Boolean to indicate whether to run single tail identification on non-normal tails, usually the dataset is so big, it is much more computationally efficient to set this as False

### Output

means: Mean of identified distributions

std: Standard deviation of identified distributions

filt: Cutoff between the distributions

## GMM.probe_filter 

GMM.probe_filter is used to filter and remove non-expressing probe sets or genes based on a predefined background threshold. The background threshold can a-priori or determined via statistical method part of this package. For example:

```
input_dataf = GMM.probe_filter(input_data,log2transform=True,filt=-0.829)
```

### Input

input_data (dataframe): input data with genes (row) x samples (columns)

log2transform (boolean): perform log2-transformation on the data

filt (float): level of background expression log2 transformed (2^-0.829) to match with the parameter log2-transform = True 

### Output 

input_dataf: return dataframe with filtered probesets 


### *Use case 2: Subcategorizing the distribution a single gene*

```
gene = 'TGFB1'

info, classif, categories,chi = GMM.GMM_modelingt(gene,input_dataf,log2transform=True,calc_backpara=True
                                    ,filt=-0.83, meanf= -3.3, stdf = 1.95)
```
### Input

gene: gene name you're interested looking at

input_dataf: input data with genes (row) x samples (columns), this can be the output of the **GMM.probe_filter** 
                    
log2transform: perform log2-transformation on the data

calc_backpara: use background threshold to draw difference between noise expression level and negative

filt: level of background expression log2 transformed 

meanf: mean of the background distribution

stdf: standard deviation of background distribution 

### Output

info: mean, covariance, and threshold of identified distribution

classif: classification annotation 

        Classification annotation:
        1: Unimodal distribution
        2: Bimodal distribution
        21: Unimodal distribution + tail
        3: Bimodal distribution + tail
        20: Biomdal distribution with chi-square fit protocol failed to fit 
        
categories: categories of each sample after GMMchisquare

chi: lowest chi-square goodness of fit value during fitting 

### *Use Case 3: Subcategorizing all genes or probes within dataset*

Below is an example of how we can use this algorithm on a large scale analysis on all genes or probes:

```
gene_name = input_dataf.index
categorize = []

for gene in tqdm(gene_name):
    info, classif, categories, chi = GMM.GMM_modelingt(gene,input_dataf, log2transform=True,calc_backpara=True
                                    ,filt=6.5924, meanf= 5.14, stdf = 1.01)
    categorize.append(categories)
    del classif, categories, chi #free up memory

    time.sleep(0.01)
```
## GMM.find_hits

This function is used to perform 2x2 contingency table analysis with the categorized data returned from the GMM.GMMmodelingt output

```
hits, filtdata = GMM.find_hits(orgi,primary='MUC2')
```
### Input

orgi: Input dataframe with categorized data that is composed of 1 or 2s (1 = low; 2 = high)

primary: Gene of interest that will be used as the primary gene compared to all other genes (index) to find correlation

### Output

Hits: 2x2 contingency table p value 

filtdata: 2x2 contingency table with p value filtered for <= 0.05

## GMM.run_hits

This function is used to output the full 2x2 contingency table from a pre-defined set of genes of interest

```
ct = GMM.run_hits(orgi,index=filtdata.index,primary='MUC2')
```

### Input

orgi: Input dataframe with categorized data that is composed of 1 or 2s (1 = low; 2 = high)

index: Pre-defined genes of interest used to compare with primary

primary: Gene of interest that will be used as the primary gene compared to all other genes (index) to find correlation 

### Output

ct: 2x2 contingency table with each component, p-value, r-value

## GMM.crosstab_table

This function is used to visualize the full 2x2 contingency table from a pre-defined set of genes of interest

```
GMM.crosstab_table(orgi,index=filtdata.index,primary='MUC2')
```
### Input

*Same as above*

# Working Example

Open python by referring to the Get Started section of this document

### Input data

| | Sample 1 | Sample 2 |
|------------ | ------------- | ------------- |
| Gene 1 | 233 | 322 |
| Gene 2 | 1022 | 333 |
| Gene 3 | 2003 | 12 |

Create a new script by clicking python3 on the drop-down menu on the top right of the shell.
Input the data in a *.csv* into python with the following:

```
import numpy as np
import pandas as np
import matplotlib.pyplot as plt
from GMMchisquare import GMMchisquare as GMM
import time
from tqdm import tqdm_notebook as tqdm

input_data = pd.read_csv(
    r'C:\Users\xxx\xxx.csv',
    index_col=[0], 
    header=0,
    na_values='---')

```

Insert your path name of your file in the C:\Users\xxx\xxx.csv. The path name of your file can be found by:

Windows: copy the file path on the top bar of the file explorer
MacOS: control + click your file, press on option to reveal copy xxx path name

### Measure the Background Threshold

You can skip this step if the a-priori suggest that there isn't a background population in your method. 
If a background population is assumed, we would want to understand whether there is a subpopulation of signals that is under the background (noise) threshold

```
means, std, filt = GMM.GMM_modelingt('Microarray Expression Data' ,input_data,log2transform=True,
                      ,calc_back=True, calc_backpara= False)
```

Now we have our mean and standard deviation of the background distribution, and the threshold seperating the two populations

### Remove Genes Below Background Threshold

If you haven't done the step above, ignore this step. 

```
input_dataf = GMM.probe_filter(input_data,log2transform=True,filt=6.5924)
```

Remember if we are log2-transforming our data here, we have to use a log2-transformed background threshold value for filt (filter threshold)

### Subcategorize Each Gene with *GMM.GMM_modelingt*

```
gene_name = input_dataf.index
categorize = []

for gene in tqdm(gene_name):
    info, classif, categories, chi = GMM.GMM_modelingt(gene,input_dataf, log2transform=True,calc_backpara=True
                                    ,filt=6.5924, meanf= 5.14, stdf = 1.01)
    categorize.append(categories)
    del classif, categories, chi #free up memory

    time.sleep(0.01)
```

This will run through every gene and categorize them using the GMMchisquare protocol. 
If there is a specific gene or probe you're interested in looking at, for example CDH1, you can change the gene_name = ['CDH1'] to look at the gene of interest or if there are multiple genes interested, for example CDH1, CDH12, CDX1, you can change the gene_name = ['CDH1', 'CDH12', 'CDX1']. 

### Saving The Output Data

```
#Create dataframe for output data
orgi = pd.DataFrame(categorize[0],index=input_dataf.columns).T

for x in tqdm(categorize[1:]):
    orgi = orgi.append(pd.DataFrame(x,index=input_dataf.columns).T)
    
    time.sleep(0.01)

#Add index for dataframe 
orgi.index = input_dataf.index

#save output at a location of choice
orgi.to_csv(r'C:/Users/xxxx/Documents/xxxx.csv')
```

Change C:/Users/xxxx/Documents/xxxx.csv to the file path you want to save your file to. 

### Analyze Categorized Data Using 2x2 Contingency Table

With all of your samples categorized, you can now run through to compare with the gene you're interested, for example MUC2. 

```
#Let's run through and filter out genes that is signficantly correlated with MUC2, a mucin marker
hits, filtdata = find_hits(orgi,primary='MUC2')

#Output contingency table using filtered index
ct = run_hits(orgi,filtdata.index,primary='MUC2')

#save output
ct.to_csv(r"C:\Users\xxxx\Documents\xxxxx.csv")
```

## Authors

* **Ta-Chun (Jeff) Liu** - [jeffliu6068](https://github.com/jeffliu6068)
* **Peter Kalugin** - *Initial work*
* **Sir Walter Fred Bodmer FRS FRSE** 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration: Thank you for all that has contributed ideas and expertise to make this possible. Let's advance science together. 
