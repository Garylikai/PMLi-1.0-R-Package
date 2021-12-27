#' PMLi: Support Functions and Dataset for Analysis of Partially Matched Samples
#'
#' Functions and dataset to support AMS 597 Project, Spring 2021, on statistical 
#' analysis strategies for partially matched samples (1st edition, 2021)
#'
#' @author Kai Li \email{kai.li@stonybrook.edu}
#'
#' @section PMLi Functions: The \code{PMLi} functions are the statistical 
#'   approaches that can be used to analyze partially matched samples.
#'
#'   \code{\link{weighted.z}} Liptak's Weighted Z-Test
#'
#'   \code{\link{modified.t}} Kim et al.'s Modified t-Statistic
#'
#'   \code{\link{corrected.z}} Looney and Jones's Corrected Z-Test
#'
#'   \code{\link{mle.hetero}}  Lin and Stivers's MLE-Based Test under
#'   Heteroscedasticity
#'
#'   \code{\link{mle.homo}} Ekbohm's MLE-Based Test under Homoscedasticity
#'   
#' @section PMLi Dataset: The \code{PMLi} sample dataset is partially matched
#'   samples.
#' 
#'   \code{\link{pm}} Sample Dataset
#'
#' @details For a complete list of functions and further details, use
#'   \code{library(help = "PMLi")}.
#'
#' @references Kuan P F, Huang B. A simple and robust method for partially
#'   matched samples using the p-values pooling approach. \emph{Statistics in
#'   medicine}. 2013; 32(19): 3247-3259.
#'
#' @docType package
#' @name PMLi
NULL
