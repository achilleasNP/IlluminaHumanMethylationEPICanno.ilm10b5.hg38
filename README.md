# IlluminaHumanMethylationEPICB5anno.ilm10b5.hg38

Illumina Human Methylation EPIC annotation 1.0 B5 using Build 38 to be used in conjunction with IlluminaHumanMethylationEPICB5manifest and minfi.

## Installation instructions
You can install it in R by using the devtools library and doing:

```r
 library(devtools)
 install_github("achilleasNP/IlluminaHumanMethylationEPICB5manifest")
 ```
and force minfi to use the package by setting the annotation of your "RGChannelSet" object as below:

## Usage instructions

```r
RGset@annotation = c(array = "IlluminaHumanMethylationEPICB5", annotation = "ilm10b5.hg38")
```
