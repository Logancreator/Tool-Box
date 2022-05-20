trim_galore -q 25 --phred33 --stringency 3 --length 20 -e 0.1 \
            --paired $dir/cmp/01raw_data/$fq1 $dir/cmp/01raw_data/$fq2  \
            --gzip -o $input_data


 # nohup htseq-count -f bam  -r name -a 25\
 # /md01/changjy/data/quail/RNAseq/RT2104103125-1/210507_SEQ053_DP8400018993TR_L01_SP2104200564/tophat_out/accepted_hits.bam\
 # /md01/changjy/data/quail/RNAseq/RT2104103126-1/210507_SEQ053_DP8400018993TR_L01_SP2104200565/tophat_out/accepted_hits.bam\
 # /md01/changjy/data/quail/RNAseq/RT2104103127-1/210507_SEQ053_DP8400018993TR_L01_SP2104200566/tophat_out/accepted_hits.bam\
 # /md01/changjy/data/quail/RNAseq/RT2104103128-1/210507_SEQ053_DP8400018993TR_L01_SP2104200567/tophat_out/accepted_hits.bam\
 # /md01/changjy/data/quail/RNAseq/RT2104103129-1/210507_SEQ053_DP8400018993TR_L01_SP2104200568/tophat_out/accepted_hits.bam\
 # /md01/changjy/data/quail/ref/ncbi_dataset/genomic.gtf\
 # >/md01/changjy/data/quail/RNAseq/count  &


featureCounts -T 15 -t exon \
-g gene_id -a ../../ref/ncbi_dataset/new_genome.gtf \
-o counts.txt LD-P10-M1 LD-P10-M2 SD-P0-M1 SD-P0-M3


library(DESeq2)
#sampleNames <- c("SD-P0-M1","SD-P0-M2","SD-P0-M3","SD-P10-M2","SD-P10-M4")
sampleNames <- c("SD-P0-M1","SD-P0-M3","SD-P10-M2","SD-P10-M4")
mydata <- read.table("featurecount.txt", header = TRUE, quote = '\t',skip =1)
names(mydata)[c(7,9,10,11)] <- sampleNames
countMatrix <- as.matrix(mydata[c(7,9,10,11)])
rownames(countMatrix) <-mydata$Geneid
#countMatrix<-countMatrix[,c(4,5,1,2,3)]
table2 <- data.frame(name = c("SD-P0-M1","SD-P0-M3","SD-P10-M2","SD-P10-M4"),condition = c("control","control","treatment","treatment"))
rownames(table2) <- sampleNames
head(countMatrix)
dds <- DESeqDataSetFromMatrix(countMatrix, colData=table2, design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds


# png("PCA.png")
# rld <- rlog(dds)
# plotPCA(rld, intgroup=c("name","condition"))
# dev.off()

# library(ggplot2)
# rld <- rlog(dds)
# data <- plotPCA(rld, intgroup=c("condition", "name"), returnData=TRUE)
# percentVar <- round(100 * attr(data, "percentVar"))
# p<- ggplot(data, aes(PC1, PC2, color=condition, shape=name)) +
# geom_point(size=3) +
# xlab(paste0("PC1: ",percentVar[1],"% variance")) +
# ylab(paste0("PC2: ",percentVar[2],"% variance"))
# p

dds <- DESeq(dds) 
results(dds)
res <- results(dds,contrast = c('condition', 'treatment', 'control'))
summary(res)
res <- res[order(res$padj),]
table(res$padj<0.05)
res <- res[order(res$padj, res$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
table(res$padj<0.05)
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
#diff_gene_deseq2 <- row.names(diff_gene_deseq2)

resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file= "/md01/changjy/data/quail/RNAseq/control_vs_treatment.csv",row.names = F)

res[which(res$log2FoldChange >= 1 & res$padj < 0.05),'sig'] <- 'up'
res[which(res$log2FoldChange <= -1 & res$padj < 0.05),'sig'] <- 'down'
res[which(abs(res$log2FoldChange) <= 1| res$padj >= 0.05),'sig'] <- 'none'
res_select <- subset(res, sig %in% c('up', 'down'))
write.table(res_select, file = "/md01/changjy/data/quail/RNAseq/control_treat.DESeq2.select.txt", sep = '\t', col.names = NA, quote = FALSE)

res_up <- subset(res, sig == 'up')
res_down <- subset(res, sig == 'down')
write.table(res_up, file = "/md01/changjy/data/quail/RNAseq/control_treat.DESeq2.up.txt", sep = '\t', col.names = NA, quote = FALSE)
write.table(res_down, file = "/md01/changjy/data/quail/RNAseq/control_treat.DESeq2.down.txt", sep = '\t', col.names = NA, quote = FALSE)

library(ggplot2)
library(ggrepel)
res<-read.csv("/md01/changjy/data/quail/RNAseq/control_vs_treatment.csv",row.names = 1)
res[which(res$log2FoldChange >= 1 & res$padj < 0.05),'sig'] <- 'up'
res[which(res$log2FoldChange <= -1 & res$padj < 0.05),'sig'] <- 'down'
res[which(abs(res$log2FoldChange) <= 1| res$padj >= 0.05),'sig'] <- 'none'

p <- ggplot(data =as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(size = 1.2) +
    scale_color_manual(values = c("#0f9b0f","gray","red"), limits = c('up', 'none', 'down')) +
    labs(x = 'log2 Fold Change', y = '-log10 adjust p-value', title = 'control vs treatment') +
    theme(plot.title = element_text(hjust = 0.8, size = 14), panel.grid = element_blank(), 
    panel.background = element_rect(color = 'black', fill = 'transparent'),

    legend.key = element_rect(fill = 'transparent')) +

    geom_vline(xintercept = c(-1, 1), lty = 4, color = 'black') + 

    geom_hline(yintercept = -log10(0.05), lty = 4, color = 'black') +

    xlim(-7, 7) + ylim(0, 45) +
    #geom_text_repel(data = res[res$padj<0.01&abs(res$log2FoldChange)>2,],aes(label = gene),size = 3,color = "black",segment.color = "black", show.legend = FALSE )
    geom_text_repel(data = res[c("TSHB","DIO2","NPY2R","Per2","Per3","DIO3","LOC107311309" ,"CGA", "GHRH", "SCL16A2"),],aes(label = gene),size = 3,color = "black",segment.color = "black", show.legend = T )
pdf('Vocalno.pdf',width = 4,height = 5)
p
dev.off()



library(org.Hs.eg.db)    
library(ggplot2)   
library(clusterProfiler) 
d1 <- read.table("./control_treat.DESeq2.down.txt", header=T, stringsAsFactor =F)
geneNames <- rownames(d1)     
#gene <-  mapIds(org.Hs.eg.db, geneNames, 'ENTREZID', 'SYMBOL')  #è¿™æ­¥è½¬æ¢éœ€è§†æƒ…å†µè€Œå®š  
BP.params <- enrichGO(   gene   = geneNames,
OrgDb  = org.Hs.eg.db,
ont   = "BP",
pAdjustMethod = "BH",
keyType = "ENSEMBL",
pvalueCutoff  = 0.01,   
qvalueCutoff  = 0.05)    
 
BP.list <- setReadable(BP.params, org.Hs.eg.db, keyType = "ENSEMBL") 



  
dotplot(BP.list, showCategory=30)
barplot(BP.list) #å¯Œé›†æŸ±å½¢å›?

 #å¯Œé›†æ°”æ³¡å›?

cnetplot(BP.list) #ç½‘ç»œå›¾å±•ç¤ºå¯Œé›†åŠŸèƒ½å’ŒåŸºå› çš„åŒ…å«å…³ç³?


heatplot(BP.list) #çƒ­å›¾å±•ç¤ºå¯Œé›†åŠŸèƒ½å’ŒåŸºå› çš„åŒ…å«å…³ç³»




goAll <- enrichGO(   gene   = geneNames, 
	OrgDb  = org.Hs.eg.db, 
	ont   = "ALL",
	pAdjustMethod = "BH",
	pvalueCutoff  = 0.01,  
	qvalueCutoff  = 0.01,
	keyType = 'ENSEMBL')

# p1 <- ggplot(data=goAll)+  
# geom_bar(aes(x=Description,y=-log10(pvalue), fill=ONTOLOGY), stat='identity') + 
# coord_flip() + 
# scale_x_discrete(limits=goAll$Description) 

# ggsave("out_bar.pdf", p1, width = 10, height=6)


# p2 <- ggplot(Edata, aes(x=GeneRatio, y=`GO description`)) +
#      geom_point(aes( size= Count , colour = -log10( pvalue ))  ) + 
#      scale_y_discrete(limits=Edata$`GO description`)+
#      ggtitle("GO enrichment")  +  
#      scale_color_gradient(low = 'green', high = 'red') + 
#      xlim(range(Edata$GeneRatio)) +
#      theme(axis.text.x=element_text(angle=0,size=8, vjust=0.7), axis.text.y=element_text(angle=0,size=6, vjust=0.7),plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16), panel.background = element_rect(fill="white", colour='gray'), panel.grid.major = element_line(size = 0.05, colour = "gray"), panel.grid.minor.y = element_line(size=0.05, colour="gray"), panel.grid.minor.x = element_line(size=0.05, colour="gray")
# )

# ggsave("out_GO.pdf", p2, width = 8, height=7)

options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")


library(DESeq2)
data=read.table("../../online_data/RNA-seq/counts.txt",header=T)
sampleNames <- c("control_1","control_2","control_3","control_4","ld28_1","ld28_2","ld28_3","ld28_4")
rownames(data)=data[,1]
countMatrix <- as.matrix(data[,c(3:10)])
table2 <- data.frame(name = c("control_1","control_2","control_3","control_4","ld28_1","ld28_2","ld28_3","ld28_4"),condition = c("control","control","control","control","treatment","treatment","treatment","treatment"))
rownames(table2) <- sampleNames
head(countMatrix)
dds <- DESeqDataSetFromMatrix(countMatrix, colData=table2, design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds
