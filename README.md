# nucleotide-diversity

## Description
This repository contains R code and data estimating the nucleotide diversity from SARS-CoV-2 wastewater samples collected in New York. 

## Citation

Dustin T. Hill, Rafael Schulman, Christopher Dunham, Yifan Zhu, Ian Vasconcellos Caldas, Yasir Ahmed-Braimah, Daryl Lamson, Lindsey Rickerman, Kirsten St. George, Hyatt Green, Brittany L. Kmush, Frank Middleton, David A. Larsen. Genetic variability of SARS-CoV-2 in wastewater and associations with community transmission. (2026). Science. Vol 392 (6799).[https://doi.org/10.1101/2025.10.24.25338735](https://doi.org/10.1101/2025.10.24.25338735](https://doi.org/10.1126/science.aed6094))

## How to use this repository
This repository contains several folders with data and R scripts. 

### *data*
This folder contains the individual and merged files used in the analysis scripts.

### *variants* 
This folder contains raw Freyja variant counts for 5 samples. These are example data provided so users can run the scripts.

### `seq diversity-1-Pi and Shannon calculation.R`
This is an R script that calculates nucleotide diversity (Pi) and Shannon diversity for individual wastewater samples that have been sequenced. The script pulls data primarily from the *variants* data.

### `seq diversity-2-Data preparation.R`
This R script loads and preprocesses the output from the `seq diversity-1-Pi and Shannon calculation.R` script. It also preprocesses concentration and clinical data that are saved in the *data* folder.

### *Figures* folder
This folder contains the R scripts for generating all the figures in the main manuscript.

### *Supplemental files* folder
This folder contains the R scripts for generating the figures and tables in the supplemental document.

### `seq diversity - functions.R`
This R script contains functions used in the analysis and supplemental scripts. You can query a help page for the functions from the `docstring` package using `?function_name`. Install `docstring` first using `install.packages("doctrinsg")`

### `seq diversity - plot themes.R`
This R script contains plot themes used in the figures.

## Equations used

### Pi equation

#### Pi per base

$\pi_{s} =( \frac{n}{n-1} ) (1 - \sum{f^{2}})$

Where $n$ is the total number of reads spanning that position, $f$ is the frequency of a variant, and the sum is over all variants at that position.


#### Windowed Pi

$\pi_{w} = \frac{1}{L} \sum{\pi_s}$

Where $L$ is the window size in base pairs. Positions with no variation in the sample are considered to have $\pi_s$ = 0. 

#### Mean Pi

The genomewide mean pi is the mean of the windowed values.

### Shannon Equation

I used the shannon equation from this link from statology: [Shannon equation](https://www.statology.org/shannon-diversity-index/)

#### Shannon per base

$H_{s} = \sum{(log(prop) * prop) * -1} $

$prop = \frac{\text{frequency of allele}}{\text{total alleles observed}} $


#### Windowed Shannon
$H_{w} = \frac{1}{L} \sum{H_{s}}$


#### Mean Shannon

The genomewide mean Shannon H is the mean of the windowed values.
