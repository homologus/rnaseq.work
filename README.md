# R library - rnaseq.work

This library provides simple APIs for RNAseq data analysis. 

## Function 'rna_diff_expr'

The function rna_diff_expr() allows analysis of RNAseq data using
various approaches (e.g. DESeq, DESeq2, edgeR, limma-voom, sleuth, bayseq,
NOIseq, EBseq, SAMseq, etc.). You just need to give the 
count_table and design_table in the form of data frames, and the rest
will be taken care of by the function.  Only DESeq2 and edgeR are implemented
so far.

~~~~~~~~~~~~
rna_diff_expr(count_table, design_table, method="DESeq")
rna_diff_expr(count_table, design_table, method="DESeq2")
rna_diff_expr(count_table, design_table, method="edgeR")
rna_diff_expr(count_table, design_table, method="limma-voom")
rna_diff_expr(count_table, design_table, method="sleuth")
rna_diff_expr(count_table, design_table, method="bayseq")
rna_diff_expr(count_table, design_table, method="NOIseq")
rna_diff_expr(count_table, design_table, method="EBseq")
rna_diff_expr(count_table, design_table, method="SAMseq")
~~~~~~~~~~~~


## Function rna_visualize

The function rna_visualize() plots RNAseq-related data in
various ways (e.g. histogram, boxplot, density plot, clustering,
MA plot, smear plot, MDS plot, PCA plot, BCV plot, volcano plot)
in both base and ggplot libraries. You just have to change
the 'method' parameter to appropriate plot type.

~~~~~~~~~~
rna_visualize(rna_data, method="hist", lib="base")
rna_visualize(rna_data, method="hist", lib="ggplot")
~~~~~~~~~~


