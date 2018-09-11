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


library(ggplot2)
library(dplyr)

#' RNAseq Visualize
#'
#' Visualize
#' https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/plotSmear
#' @export
#' @examples
#' rna_visualize(data, method="hist", lib="base")
#' rna_visualize(data, method="hist", lib="ggplot")
#' rna_visualize(data, method="boxplot", lib="base")
#' rna_visualize(data, method="boxplot", lib="ggplot")
#' rna_visualize(data, method="density", lib="base")
#' rna_visualize(data, method="density", lib="ggplot")
#' rna_visualize(data, method="cluster-gene", lib="base")
#' rna_visualize(data, method="cluster-gene", lib="ggplot")
#' rna_visualize(data, method="cluster-expt", lib="base")
#' rna_visualize(data, method="cluster-expt", lib="ggplot")
#' rna_visualize(data, method="MA", lib="base")
#' rna_visualize(data, method="MA", lib="ggplot")
#' rna_visualize(data, method="smear", lib="base")
#' rna_visualize(data, method="smear", lib="ggplot")
#' rna_visualize(data, method="MDS", lib="base")
#' rna_visualize(data, method="MDS", lib="ggplot")
#' rna_visualize(data, method="PCA", lib="base")
#' rna_visualize(data, method="PCA", lib="ggplot")
#' rna_visualize(data, method="BCV", lib="base")
#' rna_visualize(data, method="BCV", lib="ggplot")
#' rna_visualize(data, method="volcano", lib="base")
#' rna_visualize(data, method="volcano", lib="ggplot")


rna_visualize <- function(data, method="hist", lib="base"){

  if(method=="hist") {
    if(lib=="ggplot") {
      print("using hist/ggplot")
    }
    if(lib=="base") {
      print("using hist/base")
    }
  }

  if(method=="boxplot") {
    if(lib=="ggplot") {
      print("using boxplot/ggplot")
    }
    if(lib=="base") {
      print("using boxplot/base")
    }
  }

  if(method=="density") {
    if(lib=="ggplot") {
      print("using density/ggplot")
    }
    if(lib=="base") {
      print("using density/base")
    }
  }

  if(method=="cluster-expt") {
    if(lib=="ggplot") {
      print("using cluster-expt/ggplot")
    }
    if(lib=="base") {
      print("using cluster-expt/base")
    }
  }

  if(method=="cluster-gene") {
    if(lib=="ggplot") {
      print("using cluster-gene/ggplot")
    }
    if(lib=="base") {
      print("using cluster-gene/base")
    }
  }

  if(method=="smear") {
    if(lib=="ggplot") {
      print("using smear/ggplot")
    }
    if(lib=="base") {
      print("using smear/base")
      # plotSmear()
    }
  }


  if(method=="MA") {
    if(lib=="ggplot") {
      print("using MA/ggplot")
    }
    if(lib=="base") {
      print("using MA/base")
    }
  }

  if(method=="MDS") {
    if(lib=="ggplot") {
      print("using MDS/ggplot")
    }
    if(lib=="base") {
      print("using MDS/base")
      # plotMDS()
    }
  }

  if(method=="PCA") {
    if(lib=="ggplot") {
      print("using PCA/ggplot")
    }
    if(lib=="base") {
      print("using PCA/base")
      # plotPCA()
    }
  }

  if(method=="BCV") {
    if(lib=="ggplot") {
      print("using BCV/ggplot")
    }
    if(lib=="base") {
      print("using BCV/base")
      # plotBCV()
    }
  }

  if(method=="volcano") {
    if(lib=="ggplot") {
      print("using volcano/ggplot")
      # print(ggplot(data)+aes(x=log2FoldChange,y=pvalue)+geom_point())
    }
    if(lib=="base") {
      print("using volcano/base")
      # plot(data$log2FoldChange, -log10(data$pvalue), main="Volcano Plot")
    }
  }

}


