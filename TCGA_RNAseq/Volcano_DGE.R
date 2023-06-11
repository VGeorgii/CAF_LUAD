library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("survival")
library("survminer")
library("gProfileR")
library("genefilter")
library("dplyr")
library("EDASeq")
library('GO.db')
library("org.Hs.eg.db")
library ("ggplot2")
library ("EnhancedVolcano")
library ("clusterProfiler")
library ("org.Hs.eg.db")
library ("AnnotationDbi")
library ("venn")
library ("CeTF")
library ("dplyr")
library('gprofiler2')
library('Hmisc')
library('tidyr')
library('corrplot')
library('stringr')
library("PerformanceAnalytics")
library('pheatmap')
library("ComplexHeatmap")
library("RColorBrewer")
library('colorRamp2')
library("survminer")
library('survival')
library('survival')
library('data.table')



tcga_data = readRDS(file = "tcga_data.RDS")


colData(tcga_data)$tumor_stage = gsub("[ABC]$", "", colData(tcga_data)$ajcc_pathologic_stage)
colData(tcga_data)[(colData(tcga_data)[,'definition'] == 'Solid Tissue Normal'),'tumor_stage'] = 'Normal'
colData(tcga_data)[which(colData(tcga_data)$tumor_stage == "not reported"), "tumor_stage"] = NA
colData(tcga_data)$tumor_stage[is.na(colData(tcga_data)$tumor_stage)] = 'Unknown'


colData(tcga_data)[(colData(tcga_data)[,'tumor_stage'] == 'Stage I') | (colData(tcga_data)[,'tumor_stage'] == 'Stage II') | (colData(tcga_data)[,'tumor_stage'] == 'Stage III') | (colData(tcga_data)[,'tumor_stage'] == 'Stage IV'),]
#clin_df = colData(tcga_data)
clin_df = colData(tcga_data)[(colData(tcga_data)[,'tumor_stage'] == 'Stage I') | (colData(tcga_data)[,'tumor_stage'] == 'Stage II'),]

count_df = assay(tcga_data)[,rownames(clin_df)]
count_df = as.data.frame(count_df)
t_count_df = t(count_df)
clin_df[rownames(((t_count_df[t_count_df[,'ENSG00000196616.14'] > quantile(t_count_df[,'ENSG00000196616.14'])[2], ]) | (t_count_df[t_count_df[,'ENSG00000196616.14'] < quantile(t_count_df[,'ENSG00000196616.14'])[4], ]))), 'status'] = 'intermediate'
clin_df[rownames(t_count_df[t_count_df[,'ENSG00000196616.14'] >= quantile(t_count_df[,'ENSG00000196616.14'])[4], ]), 'status'] = 'high'
clin_df[rownames(t_count_df[t_count_df[,'ENSG00000196616.14'] <= quantile(t_count_df[,'ENSG00000196616.14'])[2], ]), 'status'] = 'low'



limma_pipeline = function(
    count_df,
    clin_df,
    gene_df,
    condition_variable,
    reference_group=NULL){
  
  design_factor = clin_df[, condition_variable, drop=T]
  
  group = factor(design_factor)
  if(!is.null(reference_group)){group = relevel(group, ref=reference_group)}
  
  design = model.matrix(~ group)
  
  dge = DGEList(counts=count_df,
                samples=clin_df,
                genes=gene_df)
  
  # filtering
  keep = filterByExpr(dge,design)
  dge = dge[keep,,keep.lib.sizes=FALSE]
  rm(keep)
  
  # Normalization (TMM followed by voom)
  dge = calcNormFactors(dge)
  v = voom(dge, design, plot=TRUE)
  
  # Fit model to data given design
  fit = lmFit(v, design)
  fit = eBayes(fit)
  
  # Show top genes
  topGenes = topTable(fit, coef=ncol(design), number=10000, sort.by="p")
  
  return(
    list(
      voomObj=v, # normalized data
      #fit=fit, # fitted model and statistics
      topGenes=topGenes # the 100 most differentially expressed genes
    )
  )
}

limma_res = limma_pipeline(
  count_df = count_df,
  clin_df = clin_df,
  gene_df = gene_df,
  condition_variable = "status",
  reference_group = 'low'
)


df = limma_res$topGenes


EnhancedVolcano(df, x = 'logFC', 
                y = 'adj.P.Val', 
                lab = df$gene_name,
                pointSize = 3.0,
                labSize = 4.0)