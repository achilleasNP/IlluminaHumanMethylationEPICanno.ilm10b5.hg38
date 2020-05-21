processUCSCsnp <- function(snpfile) {
    require(GenomicRanges)
    cat("Reading file\n")
    df <- read.delim(gzfile(snpfile), header = FALSE,
                     stringsAsFactors = FALSE)
    names(df) <- c("chr", "start", "end", "name", "strand",
                   "refNCBI", "class", "alleleFreqs")
    print(table(df$chr))
    cat("Only keeping chrs 1-22, X, Y\n")
    df <- df[df$chr %in% paste0("chr", c(1:22, "X", "Y")),]
    print(table(df$class))
    cat("Only keeping class 'single'\n")
    df <- df[df$class == "single",]
    cat("Computing MAF\n")
    df$alleleFreqs <- sub(",$", "", df$alleleFreqs)
    sp <- strsplit(df$alleleFreqs, ",")
    minFreq <- sapply(sp, function(xx) min(as.numeric(xx)))
    cat("Instantiating object\n")
    grSNP <- GRanges(seqnames = df$chr, strand = df$strand,
                     ranges = IRanges(start = df$start + 1, end = df$end),
                     MAF = minFreq, ref = df$refNCBI)
    names(grSNP) <- df$name
    grSNP
}

for (file in list.files("files", pattern="*_small.txt.gz")){
    objname <- paste0("gr",
		      gsub("^s", "S", unlist(strsplit(file, "_"))[[1]]),
		      "Single")
    file <- file.path("files",file)
    assign(objname, processUCSCsnp(file))
    outfile <- file.path("objects", paste(objname, "rda", sep="."))
    save(list=c(objname), file = outfile)
	
}

