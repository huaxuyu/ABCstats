
#' Adaptive Box-Cox transformation
#'
#' @description
#' Feature-specific data transformation to improve data normality in untargeted metabolomics.
#'
#' @param FeatureTable Data frame with features in row and samples in column (default).
#' @param Impt A single string specifying the imputation method to be used.
#' @param LambdaRange A numeric vector indicating the optimization range of lambda value.
#' @param GapIdentifier A single value to identify gaps in \code{FeatureTable}.
#' @param SampleInCol \code{TRUE} if samples are in column. \code{FALSE} if samples are in row.
#' @param Output \code{TRUE} will output the result table in the current working directory.
#'
#' @details
#'
#' Adaptive Box-Cox (ABC) transformation is designed to improve the data normality in untargeted metabolomics.
#' ABC transformation contains three modules: \cr
#' 1. Data pre-treatment by gap-filling and data scaling \cr
#' 2. Optimize lambda value for ABC transformation \cr
#' 3. Perform the ABC data transformation using the optimized lambda value \cr
#'
#' \code{FeatureTable} contains measured signal intensities of metabolic features,
#' with features in row and samples in column (default). The column names should
#' be sample names, and the first row should be sample group names (e.g. control, case).\cr
#'
#' Four imputation methods are provided here: \cr
#' 1. \code{default}, default imputation method by ABC transformation \cr
#' 2. \code{knn}, KNN method supported by \pkg{VIM} package. See \code{\link[VIM]{kNN}} for details. \cr
#' 3. \code{addition}, only replace gaps with identical small values (min/5). \cr
#' 4. \code{rf}, random forest method supported by \pkg{mice} package. See \code{\link[VIM]{mice}} for details.
#' This method is not recommended due to long calculation time when feature number > 200. \cr
#'
#' @return
#' This function returns the transformed data frame.
#'
#' @export
#'
#' @references To be updated.
#'
#' @examples
#' TransformedTable = ABCtransform(DemoData)


ABCtransform = function(FeatureTable, Impt="default", LambdaRange=c(-3,3),
                        GapIdentifier=0, SampleInCol=TRUE, Output=FALSE){

  message("ABC transformation is running...")

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
  # Add a new column for storing lambda value
  GroupNames = as.character(FeatureTable[1,-1])
  FeatureTable$lambda = -99
  FeatureTable$lambda[1] = "lambda"

  pb = txtProgressBar(min = 0, max = nrow(result_table), style = 3)
  # Process each individual feature
  for (feature in 1:nrow(result_table)) {
    # Extract the quantitative data
    data_seq = as.numeric(result_table[feature,])

    # ABC transformation
    opt_result = lambdaOpt(data_seq, GroupNames, L1 = LambdaRange[1], L2 = LambdaRange[2])
    FeatureTable$lambda[feature+1] = opt_result$lambda
    if(is.numeric(opt_result$lambda)){
      result_table[feature,] = opt_result$data_trans
    }
    setTxtProgressBar(pb, feature)
  }
  close(pb)

  FeatureTable[-1,-c(1, ncol(FeatureTable))] = result_table
  # Output
  if (Output) {
    write.csv(FeatureTable,"transformed_table.csv", row.names = F)
  }
  message("ABC transformation is finished.")
  return(FeatureTable)
}
