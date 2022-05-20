library(dplyr)
library(grid)
library(ggplot2)
library(gridExtra)
library(Rsamtools)
library(GenomicRanges)
options(stringsAsFactors=FALSE)
rm(list = ls())
gtf <- read.table("/md01/changjy/data/quail/ref/ncbi_dataset/quail.gtf", sep="\t")
gtf_sub <- gtf[which(gtf[, 3] == "gene"),]
#genes <- GRanges(seqnames=Rle(gtf_sub[,1]), ranges=IRanges(gtf_sub[, 4], end=gtf_sub[, 5]), strand=Rle(strand(gtf_sub[, 7])),
#gene_id=gsub("\\..*", "", gsub(".*gene_id ", "", gtf_sub[, 9])), symbol=gsub(";.*", "", gsub(".*gene_name ", "", gtf_sub[, 9])))
x=strsplit(gsub("\\..*", "", gsub(".*gene_id ", "", gtf_sub[, 9])),split=";")
genes_id=as.character(lapply(x, function(x) x[1]))
genes <- GRanges(seqnames=Rle(gtf_sub[,1]), ranges=IRanges(gtf_sub[, 4], end=gtf_sub[, 5]), strand=Rle(strand(gtf_sub[, 7])),
                 gene_id=genes_id,symbol=genes_id)

genes
gtf_sub <- gtf[which(gtf[, 3] == "exon"),]
#exons <- GRanges(seqnames=Rle(gtf_sub[,1]), ranges=IRanges(gtf_sub[, 4], end=gtf_sub[, 5]), strand=Rle(strand(gtf_sub[, 7])),
#gene_id=gsub("\\..*", "", gsub(".*gene_id ", "", gtf_sub[, 9])), symbol=gsub(";.*", "", gsub(".*gene_name ", "", gtf_sub[, 9])))

x=strsplit(gsub("\\..*", "", gsub(".*gene_id ", "", gtf_sub[, 9])),split=";")
genes_id=as.character(lapply(x, function(x) x[1]))
exons <- GRanges(seqnames=Rle(gtf_sub[,1]), ranges=IRanges(gtf_sub[, 4], end=gtf_sub[, 5]), strand=Rle(strand(gtf_sub[, 7])),
                 gene_id=genes_id,symbol=genes_id)
exons


gtf_sub <- gtf[which(gtf[, 3] == "start_codon"),]
#TSS <- GRanges(seqnames=Rle(gtf_sub[,1]), ranges=IRanges(gtf_sub[, 4], end=gtf_sub[, 4]), strand=Rle(strand(gtf_sub[, 7])),
#gene_id=gsub("\\..*", "", gsub(".*transcript_id ", "", gtf_sub[, 9])), symbol=gsub(";.*", "", gsub(".*transcript_name ", "", gtf_sub[, 9])))
x=strsplit(gsub("\\..*", "", gsub(".*gene_id ", "", gtf_sub[, 9])),split=";")
genes_id=as.character(lapply(x, function(x) x[1]))
TSS<-GRanges(seqnames=Rle(gtf_sub[,1]), ranges=IRanges(gtf_sub[, 4], end=gtf_sub[, 5]), strand=Rle(strand(gtf_sub[, 7])),
           gene_id=genes_id,symbol=genes_id)

geneAnnotation <- SimpleList(genes=genes, exons=exons, TSS=TSS)
saveRDS(geneAnnotation, "/md01/changjy/data/quail/ref/ncbi_dataset/gencode.ct22.annotation.rds")

geneAnnotation$genes
geneAnnotation$exons
geneAnnotation$TSS


#library(stringr)
#str_sub(gsub("\\..*", "", gsub(".*gene_id ", "", gtf_sub[, 9])),1,str_locate(gsub("\\..*", "", gsub(".*gene_id ", "", gtf_sub[, 9])),"\\;")[2])


get_polygon <- function(x, y, w, h, arrow="*", rate=1) {
  if (arrow == "-")
  {
    res <- data.frame(x=x, y=y)
    res <- rbind(res, c(x+min(w, rate*h/2), y+h/2))
    res <- rbind(res, c(x+w, y+h/2))
    res <- rbind(res, c(x+w, y-h/2))
    res <- rbind(res, c(x+min(w, rate*h/2), y-h/2))
    return(res)
  }
  res <- data.frame(x=x, y=y+h/2)
  if (arrow == "+") {
    res <- rbind(res, c(x+max(0, w-rate*h/2), y+h/2))
    res <- rbind(res, c(x+w, y))
    res <- rbind(res, c(x+max(0, w-rate*h/2), y-h/2))
  } else {
    res <- rbind(res, c(x+w, y+h/2))
    res <- rbind(res, c(x+w, y-h/2))
  }
  res <- rbind(res, c(x, y-h/2))
  return(res)
}
setwd("/md01/changjy/data/quail/mapping/sort.rmdup.bam/merge")

#geneAnnotation <- readRDS("/md01/changjy/data/quail/ref/ncbi_dataset/gencode.ct22.annotation.rds")
file_list <- gsub(".bam", "", list.files(path="/md01/changjy/data/quail/mapping/sort.rmdup.bam/merge", pattern="*.bam$"))
file_list

read_count <- data.frame(Sample=file_list, Count=0)
read_count
for (i in 1:length(file_list)) read_count$Count[i] <- sum(read.delim(paste0(file_list[i], "_chr.tsv"), h=F)[, 2])
read_base <- min(read_count$Count)
read_count$Rate <- floor(log10(read_count$Count/read_base))
read_count$Rate[which(read_count$Rate == 0)] <- 1
read_count
read_base

data=read.table("/md01/changjy/data/quail/RNAseq/control_treat.DESeq2.select.txt")
data
gene_list <- c("PGR")
data_info <- data.frame(ID=gene_list, Symbol=gene_list)

tileSize <- 5
expand <- 3000
imgWidth <- 12
for (geneSymbol in gene_list)
{
  print(geneSymbol)
  region <- geneAnnotation$genes[which(mcols(geneAnnotation$genes)$gene_id == geneSymbol)]
  chr <- names(table(as.character(seqnames(region))))
  strand <- names(table(as.character(strand(region))))
  if (length(chr) != 1 | length(strand) != 1) print("multiple terms")
  if(strand == "+") start(region) <- start(region)-expand else end(region) <- end(region)+expand
  pos <- c(trunc(max(min(start(region))-tileSize, 0)/tileSize)*tileSize, trunc((max(end(region))+tileSize)/tileSize)*tileSize)
  subrng <- GRanges(seqnames=chr, ranges=IRanges(pos[1], end=pos[2]))
  regionTiles <- seq(trunc(pos[1]/tileSize), trunc(pos[2]/tileSize))*tileSize
  
  frag_list <- SimpleList()
  for (file in paste0(file_list, ".bam")) frag_list <- c(frag_list, scanBam(BamViews(file, bamRanges=subrng)))
  group_mat <- matrix(0, nrow=length(regionTiles), ncol=length(file_list))
  for (i in 1:length(file_list))
  {
    ts <- match(trunc(frag_list[[i]][[1]]$pos/tileSize)*tileSize, regionTiles, nomatch = 0)
    ids <- which(ts > 0)
    te <- match(trunc((frag_list[[i]][[1]]$pos+frag_list[[i]][[1]]$isize)/tileSize)*tileSize, regionTiles, nomatch = 0)
    ide <- which(te > 0)
    for (nb in c(ts[ids], te[ide])) group_mat[nb, i] <- group_mat[nb, i] + 1
  }}

df <- data.frame(which(group_mat > 0, arr.ind=TRUE))
df$y <- group_mat[cbind(df[,1], df[,2])]
dfm1 <- df
dfm1$row <- dfm1$row - 1
dfm1$y <- 0
dfp1 <- df
dfp1$row <- dfp1$row + 1
dfp1$y <- 0
df <- rbind(df, dfm1, dfp1, data.frame(row=rep(c(1, length(regionTiles)), each=ncol(group_mat)), col=rep(1:ncol(group_mat), 2), y=0))
df <- df[!duplicated(df[,1:2]),]
df <- df[df$row > 0 & df$row < (length(regionTiles)+1),]
df$x <- regionTiles[df$row]
df$group <- file_list[df$col]
df <- df[order(df$group, df$x),]
df <- df[,c("x", "y", "group")]
df$y <- round(df$y/read_count$Rate[match(df$group, read_count$Sample)])
df$group <- factor(df$group, levels=file_list)
ylim <- c(0, quantile(df$y, probs=c(0.999)))
df$y[df$y < ylim[1]] <- ylim[1]
df$y[df$y > ylim[2]] <- ylim[2]
#write.csv(df, paste0("frag_test_", groupBy, "_", geneSymbol, ".csv"))

gene_info <- geneAnnotation$genes[which(mcols(geneAnnotation$genes)$gene_id == geneSymbol)]
mcols(gene_info)$gene_id <- paste(mcols(gene_info)$gene_id, 1:length(gene_info), sep=".")
exon_info <- data.frame()
for (i in 1:length(gene_info))
{
  exon_sub <- reduce(subsetByOverlaps(geneAnnotation$exons, gene_info[i]))
  exon_sub <- data.frame(start=start(exon_sub), end=end(exon_sub), id=mcols(gene_info)$gene_id[i], strand="*")
  if (strand == "+") exon_sub$strand[which.max(exon_sub$end)] <- "+"
  if (strand == "-") exon_sub$strand[which.min(exon_sub$start)] <- "-"
  exon_info <- rbind(exon_info, exon_sub)
}
gene_info <- data.frame(start=start(gene_info), end=end(gene_info), id=mcols(gene_info)$gene_id, strand=strand)

term <- data_info$Symbol[match(geneSymbol, data_info$ID)]
tss_min <- 0
tss_max <- 0
if (strand == "+") tss_min <- min(gene_info$start)-1000 else tss_min <- min(gene_info$end)-1000
if (strand == "+") tss_max <- max(gene_info$start)+1000 else tss_max <- max(gene_info$end)+1000
tss_min <- max(tss_min, pos[1])
tss_max <- min(tss_max, pos[2])
pa <- ggplot(df, aes(x=x, y=y))+geom_area(stat="identity", colour="gray20", fill="gray20")+facet_wrap(facets= ~ group, strip.position='left', ncol=1)+
  scale_x_continuous(expand=expansion(), limits=c(pos[1], pos[2]))+guides(fill=FALSE, colour=FALSE)+
  geom_vline(xintercept=tss_min, colour="gray60", linetype="dashed", size=1)+
  
  geom_vline(xintercept=tss_max, colour="gray60", linetype="dashed", size=1)+
  theme(panel.background=element_blank(), panel.spacing=unit(0, "lines"), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
        strip.text.y.left=element_text(size=12, color="black", angle=0, hjust=1), strip.background=element_blank())

bg <- data.frame(get_polygon(pos[1], 0, pos[2]-pos[1], nrow(gene_info)), group="Gene")
pc <- ggplot(bg, aes(x=x, y=y))+geom_polygon(color="white", fill="white", alpha=0)+facet_wrap(facets= ~ group, strip.position='left', ncol=1)+
  scale_x_continuous(expand=expansion(), limits=c(pos[1], pos[2]), breaks=c(tss_min, tss_max))+
  labs(x=paste0("\nATAC Signal Coverage (0-", max(df$y), ")\nRelated to ", term, " (", chr, ": ", pos[1], "-", pos[2], ")"))+
  geom_vline(xintercept=tss_min, colour="red", linetype="dashed", size=1)+
  #geom_vline(xintercept=6187792, colour="red", linetype="dashed", size=1)
#geom_vline(xintercept=6179270, colour="red", linetype="dashed", size=1)+
  geom_vline(xintercept=tss_max, colour="red",linetype="dashed", size=1)+
  theme(panel.background=element_blank(), panel.spacing=unit(0, "lines"), axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(),
        axis.title.x=element_text(size=15, color="black"), axis.line.x=element_line(linetype=1,colour="black"), axis.text.x=element_text(size=12, colour="black"),
        strip.text.y.left=element_text(size=12, color="black", angle=0, hjust=1), strip.background=element_blank())
for (i in 1:nrow(gene_info)) pc <- pc+geom_polygon(data=get_polygon(gene_info$start[i], i-nrow(gene_info)/2-0.5, gene_info$end[i]-gene_info$start[i], 0.1), 
                                                   aes(x=x, y=y), color="gray20", fill="gray20")
eid <-match(exon_info$id, gene_info$id)
xy_rate <- (pos[2] - pos[1])/50
for (i in 1:nrow(exon_info)) pc <- pc+geom_polygon(data=get_polygon(exon_info$start[i], eid[i]-nrow(gene_info)/2-0.5, exon_info$end[i]-exon_info$start[i], 0.8, 
                                                                    exon_info$strand[i], xy_rate), aes(x=x, y=y), color="gray20", fill="gray20")


imageHeight <- imgWidth*(length(file_list)+1+nrow(gene_info)/2)/20
ggsave(plot=patchwork::wrap_plots(A=pa, B=pc, ncol=1, heights=c(length(file_list)*2, nrow(gene_info))), 
       width=imgWidth, height=imageHeight, dpi=200, filename=paste0("bulk_test_sample_", term, "_", geneSymbol, ".png"), limitsize=F)

}

