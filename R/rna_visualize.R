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
#' @export
#' @examples
#' rna_visualize(data,method="MDSplot")
#' https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/plotSmear
#' rna_visualize(data,method="", lib="base")
#' rna_visualize(data,method="", lib="ggplot")
#'
#' rna_visualize(data,method="hist")
#' rna_visualize(data,method="boxplot")
#' rna_visualize(data,method="density")
#' rna_visualize(data,method="smear")
#' rna_visualize(data,method="cluster-expt")
#' rna_visualize(data,method="cluster-gene")
#' rna_visualize(data,method="MA")
#' rna_visualize(data,method="MDS")
#' rna_visualize(data,method="PCA")
#' rna_visualize(data,method="BCV")
#' rna_visualize(data,method="volcano")
#' rna_visualize(data,method="venn")


rna_visualize <- function(data, method="hist", lib="base"){

  if(method=="hist") {
    if(lib=="ggplot") {
    }
    if(lib=="base") {
    }
  }

  if(method=="boxplot") {
    if(lib=="ggplot") {
    }
    if(lib=="base") {
    }
  }

  if(method=="density") {
    if(lib=="ggplot") {
    }
    if(lib=="base") {
    }
  }

  if(method=="smear") {
    if(lib=="ggplot") {
    }
    if(lib=="base") {
    }
  }

  if(method=="cluster-expt") {
    if(lib=="ggplot") {
    }
    if(lib=="base") {
    }
  }

  if(method=="cluster-gene") {
    if(lib=="ggplot") {
    }
    if(lib=="base") {
    }
  }

  if(method=="MA") {
    if(lib=="ggplot") {
    }
    if(lib=="base") {
    }
  }

  if(method=="MDS") {
    if(lib=="ggplot") {
    }
    if(lib=="base") {
    }
  }

  if(method=="PCA") {
    if(lib=="ggplot") {
    }
    if(lib=="base") {
    }
  }

  if(method=="BCV") {
    if(lib=="ggplot") {
    }
    if(lib=="base") {
    }
  }

  if(method=="volcano") {
    if(lib=="ggplot") {
    }
    if(lib=="base") {
    }
  }

  if(method=="venn") {
    if(lib=="ggplot") {
    }
    if(lib=="base") {
    }
  }

}

