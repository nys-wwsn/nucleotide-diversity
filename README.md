# nucleotide-diversity

## Pi equation

### Pi per base

$\pi_{s} =( \frac{n}{n-1} ) (1 - \sum{f^{2}})$

Where $n$ is the total number of reads spanning that position, $f$ is the frequency of a variant, and the sum is over all variants at that position.


### Windowed Pi

$\pi_{w} = \frac{1}{L} \sum{\pi_s}$

Where $L$ is the window size in base pairs. Positions with no variation in the sample are considered to have $\pi_s$ = 0. 

### Mean Pi

The genomewide mean pi is the mean of the windowed values.

## Shannon Equation

### Shannon per base

$H_{s} = \sum{(log(prop) * prop) * -1} $

$prop = \frac{frequency\of\allele}{total\alleles\observed} $


### Windowed Shannon
$H_{w} = \frac{1}{L} \sum{H_{s}}$


### mean Shannon

The genomewide mean Shannon H is the mean of the windowed values.