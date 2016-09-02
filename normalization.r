library(DESeq2)
##构建DESeq2输入文件
count.origin <- read.table("read.counts.all.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE) ##读入read.counts原始文件，行为基因，列为材料，必须为un-normalized数据
lines <- as.vector(names(count.origin)[-1]) #提取列名材料列表
genes <- as.vector(count.origin[,1]) #提取行名基因列表
count <- as.matrix(count.origin[,-1]) #将read.counts数据部分转化为matrix
dimnames(count) <- list(genes,lines) #分别将基因列表和材料列表赋值给read.counts矩阵的行名和列名，构成read.counts输入矩阵(格式为matrix)
coldata <- as.data.frame(cbind(lines,"paired-end"),row.names=lines) #构建注释信息文件，注释信息包含材料名和reads种类，还可自由添加其他注释信息，行名为材料名（格式为data.frame）
names(coldata) <- c("condition","type") #为注释文件对应的列名赋值
##DESeq2处理read.counts数据
dds <- DESeqDataSetFromMatrix(countData=count,colData=coldata,design = ~ condition) #使用DESeqDataSetFromMatrix读入已构建好的read.counts矩阵和材料注释信息，design参数代表之后的差异表达分析所依据的注释信息（例：~ condition），如果该列所有元素都相同，则为~ 1
my.counts.raw <- counts(dds,normalized=FALSE)
dds <- estimateSizeFactors(dds) #使用meidian ratio方法估计输入文件的factors，只有经过此操作的数据才能进行normalization
#Construct GRangesList from gtf file
library(rtracklayer)
gff <- import("B73.AGPv3.27.gtf","gtf") ##读入B73_v3注释文件gtf
gff.flag <- mcols(gff)$gene_biotype == "protein_coding" & mcols(gff)$type == "exon" #根据gtf注释信息提取每个基因所有exon结构
gff.exon <- gff[gff.flag]
genes_list <- split(gff.exon, mcols(gff.exon)$gene_id)
#genes_list (GRangesList object) and dds (DESeqDataSet) has already been
#order by name internally. Both have the same genes and same order.
if (length(rownames(dds)) == sum(rownames(dds) == names(genes_list)))
{print("genes_list and dds have the same gene order")}
rowRanges(dds) <- genes_list
my.fpkm.from_normalized_counts <- fpkm(dds,robust = TRUE) #fpkm计算，输入文件为raw readcounts文件，robust=TRUE代表使用基于median ratio的方法进行fpkm的normalization
write.table(cbind(genes,my.fpkm.from_normalized_counts),file="normalized_counts_FPKM.txt",sep="\t",quote=FALSE,row.names=FALSE)
