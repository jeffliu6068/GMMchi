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

GMMchi.GMMmodelingt is the function that runs GMMchi on input data. Here, we define calc_back = True to specify the use of calculating background threshold.

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

filt: Cutoff between the distributions, this is the background threshold that seperates the background vs normal distribution

## Filter and Remove Non-expressing Genes

GMMchi.probe_filter is used to filter and remove non-expressing probe sets or genes based on a background threshold determined above. The background threshold can be a-priori or determined via a statistical method included in this package (shown above).

```
input_dataf = GMMchi.probe_filter(input_data_cancer, log2transform=True, filt=-0.829)
```

### Input

filt (float): the background threshold. Note that the threshold needs to match the parameter log2transform such that if log2transform = True, filt needs to be a log2transformed threshold and vice versa

### Output 

input_dataf: return dataframe with filtered probesets 

### Categorizing the Distribution a Single Gene

```
gene = 'TGFB1' #Transforming growth factor beta 1

info, classif, categories, chi, bins, f = GMMchi.GMM_modelingt(gene, input_dataf, log2transform = True,
                                            filt=-0.83, meanf= -3.3, stdf = 1.95)
```
### Input

gene: gene of interest

input_dataf: Dataframe with genes (row) x samples (columns), this is usually the postprocessed data or the output of **GMM.probe_filter** 

meanf: mean of the background distribution (retrieved from the result of **Calculating Background Threshold**)

stdf: standard deviation of background distribution (retrieved from the result of **Calculating Background Threshold**)

### Output

info: mean(s), covariance(s), and threshold(s) of the identified distribution returned as a list of list

classif: name of the category the distribution is identified as 

        Classifications:
        1) Bimodal
        2) Unimodal
        3) Categorical unimodal
        4) Unimodal with a non-normal tail
        5) Bimodal with a non-normal tail
        6) Poorly fitted bimodal
        
categories: the returned categorized data as a list

chi: the chi-square goodness of fit of the fitted model returned as a float

bins: the bins of the histogram returned as a list

f: figure of the plot returned as a matplotlib fig object

### Large-scale Categorization of the Input Data (All genes)

Below is an example of how we can use this algorithm on a large scale analysis on all genes or probe sets:

```
genes = input_dataf.index #the index of the dataframe or a list of all genes
categorize = [] #append as list of list of categorized data

for gene in tqdm(genes):
    info, classif, categories, chi, bins, f = GMMchi.GMM_modelingt(gene, input_dataf, log2transform=True,
                                                  filt=6.5924, meanf= 5.14, stdf = 1.01)
    categorize.append(categories)
    
    del classif, categories, chi #free up memory
   
categorized_df = pd.DataFrame(categorize, index = input_dataf.index, columns = input_dataf.columns)
```
## Run a 2x2 Table Analysis

GMMchi.find_hits is used to perform a 2x2 contingency table analysis with the categorized data returned from **GMM.GMMmodelingt** on the gene of interest

```
hits, significant_hits, table_sig_hits = GMMchi.find_hits(categorized_df, primary='TGFB1')
```
### Input

categorized_df: Dataframe with categorized data that is composed of 1 or 2s (1 = low; 2 = high)

primary: Gene of interest that will be used as the primary gene compared to all other genes (index) to find correlation

### Output

Hits: 2x2 contingency table p value 

significant_hits: 2x2 contingency table with p value filtered for <= 0.05

table_sig_hits: returned as a list of ['+/+','+/-','-/+','-/-','p-value','R value', 'Inclusion Criterion']

# Working Example

Please find a working example in the example folder

## Authors

* **Ta-Chun (Jeff) Liu** - [jeffliu6068](https://github.com/jeffliu6068)
* **Peter Kalugin** - *Initial work*
* **Sir Walter Fred Bodmer FRS FRSE** 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration: Thank you for all that has contributed ideas and expertise to make this possible. Let's advance science together. 
