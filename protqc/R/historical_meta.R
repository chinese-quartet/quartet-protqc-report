#' Metadata of historical data sets.
#' @format A data frame with 384 rows and 11 variables:
#' \describe{
#'   \item{omics}{dbl Default: "Shotgun-Proteomics".}
#'   \item{library}{dbl A unique ID for each sample.}
#'   \item{batch}{dbl The batch info.}
#'   \item{extraction_date}{The extraction date.}
#'   \item{instrument}{dbl The instrument for LC-MS.}
#'   \item{site}{dbl The site/lab.}
#'   \item{platform}{dbl The platform.}
#'   \item{database}{dbl The database for library searching.}
#'   \item{tool}{dbl The tool/software.}
#'   \item{sample}{dbl The sample types (D5, D6, F7, M8).}
#'   \item{replicate}{dbl The replicate number for each sample type.}
#' }
"historical_meta"