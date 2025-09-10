library(rrvgo)
library(org.Hs.eg.db)

sig_go_ids <- gse@result$ID
go_term_qval <- gse@result$qvalue

simMatrix <- calculateSimMatrix(sig_go_ids,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(go_term_qval), sig_go_ids)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

# Plotting

heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6)


scatterPlot(simMatrix, reducedTerms)


treemapPlot(reducedTerms)


wordcloudPlot(reducedTerms, min.freq=1, colors="black")
