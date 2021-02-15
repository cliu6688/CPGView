#!/usr/bin/env Rscript

library(gggenes)
library(ggplot2)
options (warn = -1)

Args <- commandArgs(trailingOnly=TRUE)
geneData = read.csv(Args[1])
subgeneData = read.csv(Args[2])
geneData$Exon <- factor(geneData$Exon, levels=c('exon1', 'exon2 (IRa)', 'exon3 (IRa)', 'exon2 (IRb)', 'exon3 (IRb)'))
cddPalette <- c("#808080","#FFFFFF","#FFFFFF","#000000","#000000")
cl <- c("#FFFFFF","#000000","#000000","#000000","#000000","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF")
pdf(Args[3], width=6.27,height=2.32)
ggplot(geneData, aes(xmin = start, xmax = end, y = strand,label = strand,fill = Exon,forward = direction))+geom_hline(yintercept = c(subgeneData$strand),color = "white",size = 1,linetype ="dotted") + geom_gene_arrow() +geom_gene_label(colour = cl,aes(label = label_label)) + facet_wrap(~ Gene, scales = "free", ncol = 1) + scale_x_continuous(breaks = c(geneData$start,geneData$end),labels =c(geneData$label_from,geneData$label_to))+ scale_y_discrete(name="\nTrans-splicing Genes\n") + theme(axis.title.y=element_text(size=8)) + theme(axis.ticks.y=element_blank())+ theme(strip.text.x = element_blank()) +  theme(plot.title = element_text(hjust = 0.5,lineheight=0.8)) + scale_fill_manual(values = cddPalette) +theme_genes() + theme(axis.text.x=element_text(angle = 16,size = 5,hjust=0.5, vjust=0.5)) + theme(legend.title = element_text(colour='white')) + guides(colour = guide_legend(order = 2),shape = guide_legend(order = 1))
dev.off()
print("---  Pdf file has been generated, please check the folder under the current file!  ---")
