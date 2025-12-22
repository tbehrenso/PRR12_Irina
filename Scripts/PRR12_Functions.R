


# -----------------------------------------------
# Add asterisk to genes that are significant

asterisk_sig_genes <- function(gene_list, result_genes, result_padj, pthreshold=0.05){
  
  padj_for_heatmap_genes <- result_padj[match(gene_list, result_genes)]
  significant_heatmap_genes <- gene_list[padj_for_heatmap_genes < pthreshold]
  
  asterisked_genes <- ifelse(gene_list %in% significant_heatmap_genes, paste0(gene_list, '*'), gene_list)
  
  return(asterisked_genes)
}

# Testing:
# testlist <- asterisk_sig_genes(rownames(df_subset), res_NEU_KO_vs_WT$X, res_NEU_KO_vs_WT$padj)
# asterisk_sig_genes(testlist, res_NEU_HET_vs_WT$X, res_NEU_HET_vs_WT$padj)



# -----------------------------------------------












