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
#' rna_visualize(rna_data, method="hist", lib="base")
#' rna_visualize(rna_data, method="hist", lib="ggplot")
#' rna_visualize(rna_data, method="boxplot", lib="base")
#' rna_visualize(rna_data, method="boxplot", lib="ggplot")
#' rna_visualize(rna_data, method="density", lib="base")
#' rna_visualize(rna_data, method="density", lib="ggplot")
#' rna_visualize(rna_data, method="cluster-gene", lib="base")
#' rna_visualize(rna_data, method="cluster-gene", lib="ggplot")
#' rna_visualize(rna_data, method="cluster-expt", lib="base")
#' rna_visualize(rna_data, method="cluster-expt", lib="ggplot")
#' rna_visualize(rna_data, method="MA", lib="base")
#' rna_visualize(rna_data, method="MA", lib="ggplot")
#' rna_visualize(rna_data, method="smear", lib="base")
#' rna_visualize(rna_data, method="smear", lib="ggplot")
#' rna_visualize(rna_data, method="MDS", lib="base")
#' rna_visualize(rna_data, method="MDS", lib="ggplot")
#' rna_visualize(rna_data, method="PCA", lib="base")
#' rna_visualize(rna_data, method="PCA", lib="ggplot")
#' rna_visualize(rna_data, method="BCV", lib="base")
#' rna_visualize(rna_data, method="BCV", lib="ggplot")
#' rna_visualize(rna_data, method="dispersion", lib="base")
#' rna_visualize(rna_data, method="dispersion", lib="ggplot")
#' rna_visualize(rna_data, method="volcano", lib="base")
#' rna_visualize(rna_data, method="volcano", lib="ggplot")


rna_visualize <- function(rna_data, method="hist", lib="base"){

  if(method=="hist") {
    if(lib=="ggplot") {
      print("using hist/ggplot")
      # Fix 'untreated1'
      print(ggplot(rna_data)+aes(x=untreated1)+geom_histogram())
    } else {
      print("using hist/base")
      hist(rna_data$untreated1)
    }
  }

  if(method=="boxplot") {
    if(lib=="ggplot") {
      print("using boxplot/ggplot")
    } else {
      print("using boxplot/base")
    }
  }

  if(method=="density") {
    if(lib=="ggplot") {
      print("using density/ggplot")
    } else {
      print("using density/base")
    }
  }

  if(method=="cluster-expt") {
    if(lib=="ggplot") {
      print("using cluster-expt/ggplot")
    } else {
      print("using cluster-expt/base")
    }
  }

  if(method=="cluster-gene") {
    if(lib=="ggplot") {
      print("using cluster-gene/ggplot")
    } else {
      print("using cluster-gene/base")
    }
  }

  if(method=="smear") {
    if(lib=="ggplot") {
      print("using smear/ggplot")
    } else {
      print("using smear/base")
      plotSmear(rna_data)
    }
  }


  if(method=="MA") {
    if(lib=="ggplot") {
      print("using MA/ggplot")
      x=rna_data[,1]
      y=rna_data[,2]
      M = x - y
      A = (x + y)/2
      df = rna_data.frame(A, M)
      ggplot(df, aes(x = A, y = M)) + geom_point()
    } else {
      print("using MA/base")
    }
  }

  if(method=="MDS") {
    if(lib=="ggplot") {
      print("using MDS/ggplot")
    } else {
      print("using MDS/base")
      library("edgeR")
      plotMDS(rna_data)
    }
  }

  if(method=="PCA") {
    if(lib=="ggplot") {
      print("using PCA/ggplot")
    } else {
      print("using PCA/base")
      library("DESeq2")
      plotPCA(rna_data)
    }
  }

  if(method=="BCV") {
    if(lib=="ggplot") {
      print("using BCV/ggplot")
    } else {
      print("using BCV/base")
      plotBCV(rna_data)
    }
  }

  if(method=="dispersion") {
    if(lib=="ggplot") {
      print("using dispersion/ggplot")
    } else {
      print("using dispersion/base")
      plotDispEsts(rna_data)
    }
  }

  if(method=="volcano") {
    if(lib=="ggplot") {
      print("using volcano/ggplot")
      print(ggplot(rna_data)+aes(x=log2FoldChange,y=-log10(pvalue))+geom_point())
    } else {
      print("using volcano/base")
      plot(rna_data$log2FoldChange, -log10(data$pvalue), main="Volcano Plot")
    }
  }

}

