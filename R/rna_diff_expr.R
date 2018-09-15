###############################################################################
# Author:		Manoj Pratim Samanta
# Last modified:	2018
#
# Copyright (C) 2018
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################



#' RNAseq Differential Expressison
#'
#' This function allows you to identify
#' differentially expressed genes.
#' @export
#' @examples
#' rna_diff_expr(count_table, design_table, method="DESeq")
#' rna_diff_expr(count_table, design_table, method="DESeq2")
#' rna_diff_expr(count_table, design_table, method="edgeR")
#' rna_diff_expr(count_table, design_table, method="limma-voom")
#' rna_diff_expr(count_table, design_table, method="sleuth")
#' rna_diff_expr(count_table, design_table, method="baySeq")
#' rna_diff_expr(count_table, design_table, method="NOISeq")
#' rna_diff_expr(count_table, design_table, method="EBSeq")
#

rna_diff_expr <- function(count_table, design_table, method="DESeq2") {

  if(method=="DESeq") {
    print("using DESeq")
    library(DESeq)
    library(DESeq2)
    y <- DESeqDataSetFromMatrix(count_table, design_table, ~condition)
    group = colData(y)$condition
    cds <- newCountDataSet(count_table, group)
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds)
    res <- nbinomTest(cds, "untreated", "treated")
  }

  if(method=="DESeq2") {
    print("using DESeq2")
    library(DESeq2)
    dds <- DESeqDataSetFromMatrix(count_table, design_table, ~condition)
    dds <- DESeq(dds)
    res <- results(dds)
  }

  if(method=="edgeR") {
    print("using edgeR")
    library(DESeq2)
    library(edgeR)
    y <- DESeqDataSetFromMatrix(count_table, design_table, ~condition)
    group = colData(y)$condition
    y = counts(y)
    dge = DGEList(counts = y, group = group)
    dge = estimateCommonDisp(dge)
    dge = estimateTagwiseDisp(dge)
    et = exactTest(dge)
    ## Extract results from edgeR analysis
    res = topTags(et)
  }

  if(method=="limma-voom") {
    print("using limma-voom")
    library(limma)
    library(DESeq2)
    y <- DESeqDataSetFromMatrix(count_table, design_table, ~condition)
    group = colData(y)$condition
    dgel2 <- DGEList(counts=count_table, group=group)
    dgel2 <- calcNormFactors(dgel2)
    design <- model.matrix(~group)
    v <- voom(dgel2,design)
    fit <- lmFit(v,design)
    fit <- eBayes(fit)
    res=topTable(fit,coef=2,n=Inf,sort="p")
  }

  if(method=="sleuth") {
    print("using sleuth")

    # so$filter_bool <- filter_bool
    # so$filter_df <- filter_df
    # so$obs_norm_filt <- data.frame()
    # so$tpm_sf <- vector()
    # so$bs_quants <- list()
    # so$bs_summary <- list()

  so <- list(
      kal = kal_list,
      kal_versions = kal_versions,
      obs_raw = obs_raw,
      sample_to_covariates = sample_to_covariates,
      bootstrap_summary = NA,
      full_formula = full_model,
      design_matrix = design_matrix,
      target_mapping = target_mapping,
      gene_mode = gene_mode,
      gene_column = aggregation_column,
      norm_fun_counts = norm_fun_counts,
      norm_fun_tpm = norm_fun_tpm,
      transform_fun_counts = transform_fun_counts,
      transform_fun_tpm = transform_fun_tpm,
      pval_aggregate = pval_aggregate
    )



    # so <- sleuth_fit(so, ~condition, 'full')
    # so <- sleuth_fit(so, ~1, 'reduced')
    # so <- sleuth_lrt(so, 'reduced', 'full')
    # sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
    # sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
    # res =sleuth_significant
    res = 0
  }

  if(method=="baySeq") {
    print("using baySeq")
    # CD <- new("countData", data = simData, replicates = replicates, groups = groups)
    # CD <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = cl)
    # CD <- getLikelihoods(CD, cl = cl, bootStraps = 3, verbose = FALSE)
    res=0
  }

  if(method=="EBSeq") {
    print("using EBSeq")
    res=0
  }

  if(method=="NOISeq") {
    print("using NOISeq")
    res=0
  }

  as.data.frame(res)
}

