# IlluminaHumanMethylationEPICB5anno.ilm10b5.hg38

Illumina Human Methylation EPIC annotation to use with IlluminaHumanMethylationEPICB5manifest and minfi.
You can install using install_github("achilleasNP/IlluminaHumanMethylationEPICB5anno.ilm10b5.hg38")
and force minfi to use the package by setting the annotation of your "RGChannelSet" as below:

RGset@annotation = c(array = "IlluminaHumanMethylationEPICB5", annotation = "ilm10b5.hg38")

