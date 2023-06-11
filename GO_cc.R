genes_to_test = rownames(df)

#genes_to_test = df_highest

genes_to_test = str_replace(genes_to_test,
                            pattern = ".[0-9]+$",
                            replacement = "")

genes_to_test

GO_results = enrichGO(gene = genes_to_test, keyType = "ENSEMBL", OrgDb = 'org.Hs.eg.db', ont = 'CC', pvalueCutoff=0.05,)
as.data.frame(GO_results)
GO_results@result = GO_results[order(-GO_results@result$Count),]
fit = plot(dotplot(GO_results))