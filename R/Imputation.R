
#' Data imputation
#'
#' @description
#' Missing value imputation by different methods.
#'
#' @param FeatureTable Data frame with features in row and samples in column (default).
#' @param Impt A single string specifying the imputation method to be used.
#' @param GapIdentifier A numeric vector indicating the optimization range of lambda value.
#' @param SampleInCol \code{TRUE} if samples are in column. \code{FALSE} if samples are in row.
#' @param Output \code{TRUE} will output the result table in the current working directory.
#'
#' @details
#' Three imputation methods are provided here: \cr
#' 1. Default imputation method by ABC transformation \cr
#' 2. KNN method supported by \pkg{VIM} package. See \code{\link[VIM]{kNN}} for details. \cr
#' 3. Replace gaps with identical small values (max/1000). \cr
#'
#' \code{FeatureTable} contains measured signal intensities of metabolic features,
#' with features in row and samples in column (default). The column names should
#' be sample names, and the first row should be sample group names (e.g. control, case).\cr
#'
#' @return
#' This function returns the imputed data frame.
#'
#' @export
#'
#' @references To be updated.
#'
#' @examples
#' ImputedTable = Imputation(DemoData)
#'
Imputation = function(FeatureTable, Impt="default", GapIdentifier=0, SampleInCol=TRUE, Output=FALSE){

  message("Imputation is running...")

  # Transpose FeatureTable if samples are in row
  if (!SampleInCol) {
    FeatureTable = t(FeatureTable)
  }

  # Change all gaps to NA, change values to be numeric
  IntTable = FeatureTable[-1,-1]
  IntTable[IntTable == GapIdentifier] = NA
  IntTable = as.data.frame(sapply(IntTable, as.numeric))

  # Generate a data frame for imputed data
  result_table = impute(IntTable, Impt)

  FeatureTable[-1,-1] = result_table
  # Output
  if (Output) {
    write.csv(FeatureTable,"imputed_table.csv", row.names = F)
  }
  message("Imputation is finished.")
  return(FeatureTable)
}
