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


library(readr)
library(dplyr)
library(edgeR)
library(limma)
#library(Glimma)


#' RNAseq Differential Expressison
#'
#' This function allows you to identify
#' differentially expressed genes.
#' @export
#' @examples
#' rna_diff_expr(count_table, design_table, method="DEseq")
#' rna_diff_expr(count_table, design_table, method="DEseq2")
#' rna_diff_expr(count_table, design_table, method="edgeR")
#' rna_diff_expr(count_table, design_table, method="limma-voom")
#' rna_diff_expr(count_table, design_table, method="sleuth")
#' rna_diff_expr(count_table, design_table, method="bayseq")
#' rna_diff_expr(count_table, design_table, method="NOIseq")
#' rna_diff_expr(count_table, design_table, method="EBseq")
#' rna_diff_expr(count_table, design_table, method="SAMseq")
#

rna_diff_expr <- function(count_table, design_table, method="DEseq2") {

  if(method=="DEseq") {
    print("using DEseq")
    # cds <- newCountDataSet(cnts, grp.idx)
    # cds <- estimateSizeFactors(cds)
    # cds <- estimateDispersions(cds)
    # deseq.res <- nbinomTest(cds, "knockdown", "control")
    # deseq.fc=deseq.res$log2FoldChange
    # names(deseq.fc)=deseq.res$id
    # sum(is.infinite(deseq.fc))
    # deseq.fc[deseq.fc>10]=10
    # deseq.fc[deseq.fc< -10]=-10
    # exp.fc=deseq.fc
    # out.suffix="deseq"
  }

  if(method=="DEseq2") {
    library(DESeq2)
    print("using DEseq2")
    dds <- DESeqDataSetFromMatrix(count_table, design_table, ~condition)
    dds <- DESeq(dds)
    res <- results(dds)
  }

  if(method=="edgeR") {
    print("using edgeR")
    # y <- DGEList(count_table)
    # y <- calcNormFactors(y)
    # y <- estimateDisp(y, design_table)
    # fit <- glmFit(y, design_table)
    # res <- glmLRT(fit, coef = 2)

  }

  if(method=="limma-voom") {
    print("using limma-voom")
    # v <- voom(count_table, design_table, plot = TRUE)
    # fit <- lmFit(v)
    # cont.matrix <- makeContrasts(B.PregVsLac=basal.pregnant - basal.lactate,levels=design)
    # fit.cont <- contrasts.fit(fit, cont.matrix)
    # fit.cont <- eBayes(fit.cont)
    # summa.fit <- decideTests(fit.cont)
  }

  if(method=="sleuth") {
    print("using sleuth")
  }

  if(method=="bayseq") {
    print("using bayseq")
    # CD <- new("countData", data = simData, replicates = replicates, groups = groups)
    # CD <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = cl)
    # CD <- getLikelihoods(CD, cl = cl, bootStraps = 3, verbose = FALSE)
  }

  if(method=="EBseq") {
    print("using EBseq")
  }

  if(method=="NOIseq") {
    print("using NOIseq")
  }

  if(method=="SAMseq") {
    print("using SAMseq")
  }

  as.data.frame(res)
}
