# The code for the manifest creation is based on the code included in the package
# IlluminaHumanMethylationEPICanno.ilm10b4.hg19.
# The raw data used are the input files listed below with the md5sums.:
# de6945904b5b1d750ff5b76dba0b0840  MethylationEPIC_v-1-0_B5.csv
# 3bb0678989318410489ce77173c7d236  minfiDataEPIC/inst/extdata/200144450021/200144450021_R05C01_Grn.idat
# The idat file is from the minfiDataEpic package while the manifest file is from illumina 
# (https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html)
library(minfi)
manifest.filepath <- "/restricted/projectnb/fhs-methylation/resources/manifest/MethylationEPIC_v-1-0_B5.csv"

if(!file.exists(manifest.filepath) || !file.exists("objects")) {
    cat("Missing files, quitting\n")
    q(save = "no")
}

maniTmp <- minfi:::read.manifest.EPIC(manifest.filepath)
anno <- maniTmp$manifest
manifestList <- maniTmp$manifestList

## Checking
library(illuminaio)
idat.filepath = "../../../minfiDataEPIC/inst/extdata/200144450018/200144450018_R04C01_Grn.idat"
epic <- readIDAT(idat.filepath)
address.epic <- as.character(epic$MidBlock)
dropCpGs <- anno$Name[anno$AddressB != "" & !anno$AddressB %in% address.epic]
dropCpGs <- anno$Name[anno$AddressA != "" & !anno$AddressA %in% address.epic]
table(substr(dropCpGs, 1,2))


## Manifest package
IlluminaHumanMethylationEPICB5manifest <- do.call(IlluminaMethylationManifest,
                                                list(TypeI = manifestList$TypeI,
                                                     TypeII = manifestList$TypeII,
                                                     TypeControl = manifestList$TypeControl,
                                                     TypeSnpI = manifestList$TypeSnpI,
                                                     TypeSnpII = manifestList$TypeSnpII,
                                                     annotation = "IlluminaHumanMethylationEPICB5"))
## Annotation package
anno$IlmnID <- NULL
nam <- names(anno)
names(nam) <- nam
nam[c("AddressA_ID", "AddressB_ID", "AlleleA_ProbeSeq", "AlleleB_ProbeSeq",
            "Infinium_Design_Type", "Next_Base", "Color_Channel")] <-  c("AddressA", "AddressB",
                                                                         "ProbeSeqA", "ProbeSeqB",
                                                                         "Type", "NextBase", "Color")

names(nam) <- NULL
names(anno) <- nam
rownames(anno) <- anno$Name
anno <- anno[getManifestInfo(IlluminaHumanMethylationEPICB5manifest, type = "locusNames"),]

Locations <- anno[, c("CHR", "MAPINFO")]
names(Locations) <- c("chr", "pos")
Locations$pos <- as.integer(Locations$pos)
Locations$chr <- paste("chr", Locations$chr, sep = "")
Locations$strand <- ifelse(anno$Strand == "F", "+", "-")
table(Locations$chr, exclude = NULL)
rownames(Locations) <- anno$Name
Locations <- as(Locations, "DataFrame")

Manifest <- anno[, c("Name", "AddressA", "AddressB",
                     "ProbeSeqA", "ProbeSeqB", "Type", "NextBase", "Color")]
Manifest <- as(Manifest, "DataFrame")

Islands.UCSC <- anno[, c("UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island")]
names(Islands.UCSC) <- c("Islands_Name", "Relation_to_Island")
Islands.UCSC <- as(Islands.UCSC, "DataFrame")
Islands.UCSC$Relation_to_Island[Islands.UCSC$Relation_to_Island == ""] <- "OpenSea"
table(Islands.UCSC$Relation_to_Island, exclude = NULL)

SNPs.Illumina <- anno[, c("SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")]
SNPs.Illumina <- as(SNPs.Illumina, "DataFrame")

usedColumns <- c(names(Manifest), names(SNPs.Illumina), 
                 c("CHR", "MAPINFO", "Strand",
                   "Chromosome_36", "Coordinate_36", "Genome_Build"),
                 c("UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island"))
Other <- anno[, setdiff(names(anno), usedColumns)]
nam <- names(Other)
nam <- sub("_NAME", "_Name", nam)
nam[nam == "X450k_Enhancer"] <- "Methyl450_Enhancer"
nam
Other <- as(Other, "DataFrame")

## We now use an exisitng grSnp object containing a GRanges of relevant SNPs.
## This is created in a separate script

##
## SNP overlap
##

map <- cbind(Locations, Manifest)
map <- GRanges(seqnames = map$chr, ranges = IRanges(start = map$pos, width = 1),
               Strand = map$strand, Type = map$Type)
map <- minfi:::.getProbePositionsDetailed(map)
names(map) <- rownames(Locations)

## dbSNP
snp.objects <- c()
for (file in list.files("objects", pattern="*Single.rda")){
    full.path <- file.path("objects", file)
    load(full.path)
    original.objname <- gsub("\\.rda", "", file)
    objname <- gsub("grSnp", "SNPs.", original.objname)
    snp.objects <- c(snp.objects, objname)
    assign(objname, 
	   minfi:::.doSnpOverlap(map, get(original.objname)))
}



annoNames <- c("Locations", "Manifest", "SNPs.Illumina", "Islands.UCSC", "Other",
snp.objects)

for(nam in annoNames) {
    cat(nam, "\n")
    save(list = nam, file = file.path("../../data", paste(nam, "rda", sep = ".")), compress = "xz")
}
annoStr <- c(array = "IlluminaHumanMethylationEPICB5",
             annotation = "ilm10b5",
             genomeBuild = "hg38")
defaults <- c("Locations", "Manifest", "SNPs.141CommonSingle", "Islands.UCSC", "Other")
pkgName <- sprintf("%sanno.%s.%s", annoStr["array"], annoStr["annotation"],
                    annoStr["genomeBuild"])

annoObj <- IlluminaMethylationAnnotation(objectNames = annoNames, annotation = annoStr,
                              defaults = defaults, packageName = pkgName)

assign(pkgName, annoObj)
save(list = pkgName,
     file = file.path("../../data", paste(pkgName, "rda", sep = ".")), compress = "xz")
sessionInfo()
q(save = "no")




