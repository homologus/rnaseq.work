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
#'

rna_diff_expr <- function(count_table, design_table, method="DESeq2") {

  library(dplyr)

  ############################
  # DESeq implementation
  ############################
  if(method=="DESeq") {
    print("using DESeq")
    library(DESeq)
    group=design_table$condition
    cds <- newCountDataSet(count_table, group)
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds)
    levels <- unique(group)
    res <- nbinomTest(cds, levels[1], levels[2])
    res <- as.data.frame(res) %>% dplyr::mutate(FinalP=padj,logFC=log2FoldChange) %>% dplyr::select(id,logFC,FinalP)
  }

  ############################
  # DESeq2 implementation
  ############################
  if(method=="DESeq2") {
    print("using DESeq2")
    library(DESeq2)
    dds <- DESeqDataSetFromMatrix(count_table, design_table, ~condition)
    dds <- DESeq(dds)
    res <- results(dds)
    res <- as.data.frame(res) %>% tibble::rownames_to_column() %>% dplyr::mutate(id=rowname, FinalP=padj,logFC=log2FoldChange) %>% dplyr::select(id,logFC,FinalP)
  }

  ############################
  # edgeR implementation
  ############################
  if(method=="edgeR") {
    print("using edgeR")
    library(edgeR)
    group=design_table$condition
    dge = DGEList(counts = count_table, group = group)
    dge = estimateCommonDisp(dge)
    dge = estimateTagwiseDisp(dge)
    et = exactTest(dge)
    ## Extract results from edgeR analysis
    res = et$table
    res <- as.data.frame(res) %>% tibble::rownames_to_column() %>% dplyr::mutate(id=rowname, FinalP=PValue) %>% dplyr::select(id,logFC,FinalP)
  }

  ############################
  # limma-voom implementation
  ############################
  if(method=="limma-voom") {
    print("using limma-voom")
    library(limma)
    library(edgeR)
    group=design_table$condition
    dgel2 <- DGEList(counts=count_table, group=group)
    dgel2 <- calcNormFactors(dgel2)
    design <- model.matrix(~group)
    v <- voom(dgel2,design)
    fit <- lmFit(v,design)
    fit <- eBayes(fit)
    res=topTable(fit,coef=2,n=Inf)
    res <- as.data.frame(res) %>% tibble::rownames_to_column() %>% dplyr::arrange(rowname) %>% dplyr::mutate(id=rowname, FinalP=adj.P.Val) %>% dplyr::select(id,logFC,FinalP)
  }

  ############################
  # sleuth (not implemented)
  ############################
  if(method=="sleuth") {
    print("using sleuth")

    # so$filter_bool <- filter_bool
    # so$filter_df <- filter_df
    # so$obs_norm_filt <- data.frame()
    # so$tpm_sf <- vector()
    # so$bs_quants <- list()
    # so$bs_summary <- list()

    # so <- sleuth_fit(so, ~condition, 'full')
    # so <- sleuth_fit(so, ~1, 'reduced')
    # so <- sleuth_lrt(so, 'reduced', 'full')
    # sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
    # sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
    # res =sleuth_significant
    res = 0
  }

  ############################
  # baySeq (not implemented)
  ############################
  if(method=="baySeq") {
    print("using baySeq")
    #library(baySeq)
    #data(simData)
    #simData[1:10,]
    #cl <- NULL
    # condition=c("untreated","untreated","untreated","untreated","treated","treated","treated")
    #replicates <- c("simA", "simA", "simA", "simA", "simA", "simB", "simB", "simB", "simB", "simB")
    #groups <- list(NDE = c(1,1,1,1,1,1,1,1,1,1), DE = c(1,1,1,1,1,2,2,2,2,2))
    #CD <- new("countData", data = simData, replicates = replicates, groups = groups)
    #libsizes(CD) <- getLibsizes(CD)
    #CD@annotation <- data.frame(name = paste("count", 1:1000, sep = "_"))
    #CD <- getPriors.NB(CD, samplesize = 1000, estimation = "QL", cl = cl)
    #CD <- getLikelihoods(CD, cl = cl, bootStraps = 3, verbose = FALSE)
    #res=topCounts(CD, group = "DE")
    res=0
  }

  ############################
  # EBSeq (not implemented)
  ############################
  if(method=="EBSeq") {
    print("using EBSeq")
    res=0
  }

  ############################
  # NOISeq (not implemented)
  ############################
  if(method=="NOISeq") {
    print("using NOISeq")
    res=0
  }

  #as.data.frame(res)
  res
}


