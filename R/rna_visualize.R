##############################################################################
# Author:    Manoj Pratim Samanta
# Last modified:  2018
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
##############################################################################


#' RNAseq Visualize
#'
#' Visualize
#' https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/plotSmear
#' @export
#' @examples
#' rna_visualize(count_table, method="hist")
#' rna_visualize(count_table, method="counts")
#' rna_visualize(count_table, method="MA")
#' rna_visualize(count_table, design_table=design_table, method="smear")
#' rna_visualize(count_table, design_table=design_table, method="MDS")
#' rna_visualize(count_table, design_table=design_table, method="PCA")
#' rna_visualize(count_table, design_table=design_table, method="BCV")
#' rna_visualize(count_table, design_table=design_table, method="dispersions")
#' rna_visualize(count_table, design_table=design_table, method="sparsity")
#' rna_visualize(data_table, method="volcano")
#'
#' rna_visualize(data_table, method="hist", lib="base")
#' rna_visualize(data_table, method="hist", lib="ggplot")


rna_visualize <- function(data_table, method="hist", lib="base", col, design_table, gene){

  library("ggplot2")

  #######################
  # histogram
  #######################

  if(method=="hist") {
    if(lib=="ggplot" || lib=="ggplot2") {
      print("using hist/ggplot")
      pseudocount=log2(data_table+1)
      print(ggplot(pseudocount)+aes_string(x=col)+geom_histogram())
    } else {
      print("using hist/base")
      hist(log2(data_table[,col]+1))
    }
  }

  #######################
  # counts
  #######################

  if(method=="counts") {
    library("DESeq2")
    if(lib=="ggplot" || lib=="ggplot2") {
      print("using hist/ggplot")
      pseudocount=log2(data_table+1)
      print(ggplot(pseudocount)+aes_string(x=col)+geom_histogram())
    } else {
      print("using hist/base")
      dds <- DESeqDataSetFromMatrix(count_table, design_table, ~condition)
      plotCounts(dds,gene)
    }
  }


  #######################
  # smear
  #######################
  if(method=="smear") {
    if(lib=="ggplot" || lib=="ggplot2") {
      print("using smear/ggplot")
    } else {
      print("using smear/base")
      library("edgeR")
      d <- DGEList(data_table, group = design_table$condition)
      d <- estimateCommonDisp(d)
      plotSmear(d)
    }
  }


  #######################
  # MA
  #######################
  if(method=="MA") {
    if(lib=="ggplot" || lib=="ggplot2") {
      print("using MA/ggplot")
      x=log2(data_table[,col[1]] + 1)
      y=log2(data_table[,col[2]] + 1)
      M = y - x
      A = (x + y)/2
      df = data.frame(A, M)
      print(ggplot(df, aes(x = A, y = M)) + geom_point())
    } else {
      print("using MA/base")
      library("limma")
      maPlot(data_table[,col[1]], data_table[,col[2]])
    }
  }


  #######################
  # MDS
  #######################
  if(method=="MDS") {
    library("PoiClaClu")
    if(lib=="ggplot" || lib=="ggplot2") {
      #print("using MDS/ggplot")
      #poisd <- PoissonDistance(t(counts(dds)))
      #samplePDM <- as.matrix( poisd$dd )
      #rownames(samplePDM) <- paste( dds$dex, dds$cell, sep=" - " )
      #colnames(samplePDM) <- NULL
      #mds <- as.data.frame(colData(dds)) %>% cbind(cmdscale(samplePDM))
      #print(ggplot(mds, aes(x = `1`, y = `2`)) + geom_point(size = 3) + coord_fixed())
    } else {
      print("using MDS/base")
      library("edgeR")
      d <- DGEList(data_table, group = design_table$condition)
      d <- calcNormFactors(d)
      plotMDS(d,labels=design_table$condition)
    }
  }


  #######################
  # PCA
  #######################
  if(method=="PCA") {
    library("DESeq2")
    if(lib=="ggplot" || lib=="ggplot2") {
      print("using PCA/ggplot")
    } else {
      print("using PCA/base")
      dds <- DESeqDataSetFromMatrix(count_table, design_table, ~condition)
      rld=rlog(dds)
      plotPCA(rld)
    }
  }


  #######################
  # BCV
  #######################
  if(method=="BCV") {
    if(lib=="ggplot" || lib=="ggplot2") {
      print("using BCV/ggplot")
    } else {
      library("edgeR")
      print("using BCV/base")
      d <- DGEList(data_table, group = design_table$condition)
      d <- estimateCommonDisp(d)
      d <- estimateTrendedDisp(d)
      d <- estimateTagwiseDisp(d)
      plotBCV(d)
    }
  }


  #######################
  # dispersion
  #######################
  if(method=="dispersion") {
    if(lib=="ggplot" || lib=="ggplot2") {
      print("using dispersion/ggplot")
    } else {
      print("using dispersion/base")
      dds <- DESeqDataSetFromMatrix(count_table, design_table, ~condition)
      dds <- makeExampleDESeqDataSet()
      dds <- estimateSizeFactors(dds)
      dds <- estimateDispersions(dds)
      plotDispEsts(dds)
    }
  }


  #######################
  # volcano
  #######################
  if(method=="volcano") {
    if(lib=="ggplot" || lib=="ggplot2") {
      print("using volcano/ggplot")
      print(ggplot(data_table)+aes(x=log2FoldChange,y=-log10(pvalue))+geom_point())
    } else {
      print("using volcano/base")
      plot(data_table$log2FoldChange, -log10(data$pvalue), main="Volcano Plot")
    }
  }

}


