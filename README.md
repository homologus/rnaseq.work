# R library - rnaseq.work

This library provides simple APIs for RNAseq data analysis. 

## Function 'rna_diff_expr'

The function rna_diff_expr() allows analysis of RNAseq data using
various approaches (e.g. DESeq, DESeq2, edgeR, limma-voom, sleuth, baySeq,
NOISeq, EBSeq, etc.). You just need to give the 
count_table and design_table in the form of data frames, and the rest
will be taken care of by the function.  Only DESeq2 and edgeR are implemented
so far.

~~~~~~~~~~~~
rna_diff_expr(count_table, design_table, method="DESeq")
rna_diff_expr(count_table, design_table, method="DESeq2")
rna_diff_expr(count_table, design_table, method="edgeR")
rna_diff_expr(count_table, design_table, method="limma-voom")
rna_diff_expr(count_table, design_table, method="sleuth")
rna_diff_expr(count_table, design_table, method="baySeq")
rna_diff_expr(count_table, design_table, method="NOISeq")
rna_diff_expr(count_table, design_table, method="EBSeq")
~~~~~~~~~~~~


## Function rna_visualize

The function rna_visualize() plots RNAseq-related data in
various ways (e.g. histogram, boxplot, density plot, clustering,
MA plot, smear plot, MDS plot, PCA plot, BCV plot, volcano plot)
in both base and ggplot libraries. You just have to change
the 'method' parameter to appropriate plot type.

~~~~~~~~~~
rna_visualize(data_table, method="hist", col="treated1")
rna_visualize(data_table, method="MA", col=c("treated1", "untreated1"))
rna_visualize(data_table, method="counts", gene="gene1")
rna_visualize(data_table, method="PCA")
rna_visualize(data_table, method="MDS")
rna_visualize(data_table, method="BCV")
rna_visualize(data_table, method="smear")
rna_visualize(data_table, method="dispersions")
rna_visualize(data_table, method="sparsity")
rna_visualize(data_table, method="hist", lib="base", options="treated1")
rna_visualize(data_table, method="hist", lib="ggplot", options="untreated1")
~~~~~~~~~~


