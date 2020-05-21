library(minfi)
manifestFile <- "../../../data_files/from_website_02242016/MethylationEPIC_v-1-0_B1.csv"
maniTmp <- minfi:::read.manifest.EPIC(manifestFile)

## manifestFile <- "../../../data_files/from_website_11122015/MethylationEPIC_15073387_v-1-0.csv"
## maniTmpNew <- minfi:::read.manifest.EPIC(manifestFile)

## ## Comparison
## intersect(names(maniTmp$manifest), names(maniTmpNew$manifest))
## setdiff(names(maniTmp$manifest), names(maniTmpNew$manifest))
## setdiff(names(maniTmpNew$manifest), names(maniTmp$manifest))

## anno <- maniTmp$manifest
## rownames(anno) <- anno$Name
## anno <- anno[,-15]
## annoNew <- maniTmpNew$manifest
## rownames(annoNew) <- annoNew$Name
## common <- intersect(anno$Name, annoNew$Name)
## anno <- anno[common,]
## annoNew <- annoNew[common,]
## equal <- sapply(1:ncol(anno), function(ii) all(anno[,ii] == annoNew[,ii]))
## notEqualNames <- names(anno)[which(!equal)]
## tmp = sapply(notEqualNames, function(nam) sum(anno[, nam] != annoNew[, nam]))
## nam <- notEqualNames[6]
## nam
## wh <- which(anno[, nam] != annoNew[, nam])
## cbind(Old = anno[head(wh), nam], New = annoNew[head(wh), nam])

## cbind(anno[head(wh), notEqualNames[c(4,5)]], annoNew[head(wh), notEqualNames[c(4,5)]], anno[head(wh), c("CHR", "MAPINFO")])
## cbind(anno[head(wh), notEqualNames[c(4,5)]], annoNew[head(wh), notEqualNames[c(4,5)]], anno[head(wh), c("CHR", "MAPINFO")])


## file.exists(manifestFile)

maniTmp <- minfi:::read.manifest.EPIC(manifestFile)
anno <- maniTmp$manifest
manifestList <- maniTmp$manifestList

## Checking
library(illuminaio)
epic <- readIDAT("../../../data_files/Demo_Data_EPIC/200144450018/200144450018_R04C01_Grn.idat")
address.epic <- as.character(epic$MidBlock)
dropCpGs <- anno$Name[anno$AddressB != "" & !anno$AddressB %in% address.epic]
dropCpGs <- anno$Name[anno$AddressA != "" & !anno$AddressA %in% address.epic]
table(substr(dropCpGs, 1,2))


## Manifest package
IlluminaHumanMethylationEPICmanifest <- do.call(IlluminaMethylationManifest,
                                                list(TypeI = manifestList$TypeI,
                                                     TypeII = manifestList$TypeII,
                                                     TypeControl = manifestList$TypeControl,
                                                     TypeSnpI = manifestList$TypeSnpI,
                                                     TypeSnpII = manifestList$TypeSnpII,
                                                     annotation = "IlluminaHumanMethylationEPIC"))
save(IlluminaHumanMethylationEPICmanifest, compress = "xz",
     file = "IlluminaHumanMethylationEPICmanifest.rda")

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
anno <- anno[getManifestInfo(IlluminaHumanMethylationEPICmanifest, type = "locusNames"),]

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
map <- minfi:::getProbePositionsDetailed(map)
names(map) <- rownames(Locations)

## dbSNP
load("../../../snps/objects/grSnp144CommonSingle.rda")
SNPs.144CommonSingle <- minfi:::.doSnpOverlap(map, grSnp144CommonSingle)
load("../../../snps/objects/grSnp142CommonSingle.rda")
SNPs.142CommonSingle <- minfi:::.doSnpOverlap(map, grSnp142CommonSingle)
load("../../../snps/objects/grSnp141CommonSingle.rda")
SNPs.141CommonSingle <- minfi:::.doSnpOverlap(map, grSnp141CommonSingle)
load("../../../snps/objects/grSnp138CommonSingle.rda")
SNPs.138CommonSingle <- minfi:::.doSnpOverlap(map, grSnp138CommonSingle)
load("../../../snps/objects/grSnp137CommonSingle.rda")
SNPs.137CommonSingle <- minfi:::.doSnpOverlap(map, grSnp137CommonSingle)
load("../../../snps/objects/grSnp135CommonSingle.rda")
SNPs.135CommonSingle <- minfi:::.doSnpOverlap(map, grSnp135CommonSingle)
load("../../../snps/objects/grSnp132CommonSingle.rda")
SNPs.132CommonSingle <- minfi:::.doSnpOverlap(map, grSnp132CommonSingle)

annoStr <- c(array = "IlluminaHumanMethylationEPIC",
             annotation = "ilmn10b",
             genomeBuild = "hg19")
defaults <- c("Locations", "Manifest",
              "SNPs.144CommonSingle", 
              "Islands.UCSC", "Other")
annoObj <-
    IlluminaMethylationAnnotation(list(Locations = Locations,
                                       Manifest = Manifest,
                                       Islands.UCSC = Islands.UCSC,
                                       Other = Other,
                                       SNPs.Illumina = SNPs.Illumina,
                                       SNPs.144CommonSingle = SNPs.144CommonSingle,
                                       SNPs.142CommonSingle = SNPs.142CommonSingle,
                                       SNPs.141CommonSingle = SNPs.141CommonSingle,
                                       SNPs.138CommonSingle = SNPs.138CommonSingle,
                                       SNPs.137CommonSingle = SNPs.137CommonSingle,
                                       SNPs.135CommonSingle = SNPs.135CommonSingle,
                                       SNPs.132CommonSingle = SNPs.132CommonSingle
                                       ),
                                  annotation = annoStr, defaults = defaults)
validObject(annoObj)

                                       

annoName <- sprintf("%sanno.%s.%s", annoStr["array"], annoStr["annotation"],
                    annoStr["genomeBuild"])
cat("creating object:", annoName, "\n")
assign(annoName, annoObj)
save(list = annoName,
     file = paste(annoName, "rda", sep = "."), compress = "xz")

