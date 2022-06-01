library(GenomicFeatures)
library(ChIPseeker)
library(ggplot2)
library(dplyr)
library(getopt)
BiocManager::install(c("GenomicFeatures","ChIPseeker","getopt"))

spec <- matrix(
  c("gtf",  "g", 2, "character", "gtf file needed",
    "peakFile", "p", 2, "character",  "peak file needed",
    "outpath",  "o", 2, "character",  "outpath",
    "help",   "h", 0, "logical",  "Help",
    "width",   "w", 1, "logical",  "width to compute"),
  byrow=TRUE, ncol=5)
spec

# 使用getopt方法
opt <- getopt(spec=spec)
if( !is.null(opt$help) || is.null(opt$gtf) || is.null(opt$peakFile) ){
    # ... 这里你也可以自定义一些东放在里面
    cat(print(spec))
    quit()
}

if( is.null(opt$width)){
    width=opt$width
} else {
    width=2000
}



#Annotation
gtf=opt$gtf
peakFile=opt$peakFile
outpath=opt$outpath
width=2000
gtf="C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.gtf"
peakFile="C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\DAP-seq\\NR51_combined_peaks.narrowPeak"
out_path="C:\\Users\\Jeff\\OneDrive\\桌面\\Goose_ATAC_info\\DAP-seq"
gtf_path=read.table(gtf,header=F,sep="\t")
txdb <- makeTxDbFromGFF(gtf)
nanog <- readPeakFile(peakFile)
peakAnno <- annotatePeak(nanog, 
                         tssRegion=c(-width, width),
                         TxDb=txdb)
peakAnno
#用/把字符串分开，fixed为是否使用正则表达式 ```
topics = strsplit(peakFile, ".", fixed= T)[1]
topics
pdf(paste(out_path,"\\",topics,".pdf",sep=""))
#pie图
plotAnnoPie(peakAnno)
vennpie(peakAnno)
#TSS的强度分布图
peakHeatmap(f, weightCol="V5", TxDb=txdb,
            upstream=width, downstream=width,
            color=rainbow(length(nanog)))
dev.off()