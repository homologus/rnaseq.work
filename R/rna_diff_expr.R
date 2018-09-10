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
library(DESeq2)
library(edgeR)
library(limma)
#library(Glimma)


#' RNAseq Differential Expressison
#'
#' This function allows you to identify
#' differentially expressed genes.
#' @export
#' @examples
#' rna_diff_exp(method="DEseq")
#' rna_diff_exp(method="DEseq2")
#' rna_diff_exp(method="edgeR")
#' rna_diff_exp(method="limma-voom")
#' rna_diff_exp(method="sleuth")
#' rna_diff_exp(method="NOIseq")
#' rna_diff_exp(method="bayseq")
#' rna_diff_exp(method="EBseq")
#' rna_diff_exp(method="SAMseq")
#

rna_diff_expr <- function(count_file, design_file, method="DEseq") {

    if(method=="DESeq") {

    }

    if(method=="DEseq2") {
    }

    if(method=="edgeR") {
    }

    if(method=="limma-voom") {
    }

    if(method=="sleuth") {

    }

    if(method=="bayseq") {

    }

    if(method=="EBseq") {

    }

    if(method=="NOIseq") {

    }

    if(method=="SAMseq") {

    }

}


