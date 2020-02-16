library(flowCore)
library(foreach)
#' Pool a flowSet of multiple flowFrames into a single pooled flowFrame.
#'
#' @description
#' `pool_fs()` lets you pool cells from multiple flowFrames contained in
#' a `flowSet` into a single `flowFrame`. If `copy_desc = TRUE` it also
#' attempts to copy the keyword values that are common in all constituent
#' flowFrames.
#'
#' @importFrom flowCore flowFrame
#' @importFrom flowCore flowSet
#' @importFrom flowCore exprs
#' @importFrom flowCore parameters
#' @importFrom flowCore description
#' @importFrom foreach foreach
#' @importFrom foreach %do%
#' @param fs Input flowSet.
#' @param copy_desc (Optional): if `TRUE` will attempt to copy common keywords
#' to output `flowFrame`.
#' @return A pooled flowFrame.
#' @note The input fs can either be a `flowCore::flowSet` or `ncdfFlow::ncdfFlowSet`.
#' @examples
#'# Mock up two flowFrames
#' set.seed(1234)
#' library(flowCore)
#' rstrings <- function(n = 10) {
#'   a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
#'   paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
#' }
#' ff1_exprs <- matrix(rnorm(1:10000), nrow = 1000)
#' colnames(ff1_exprs) <- rstrings(dim(ff1_exprs)[[2]])
#' ff1 <- flowFrame(ff1_exprs)
#' ff2_exprs <- matrix(rnorm(1:10000), nrow = 1000)
#' colnames(ff2_exprs) <- colnames(ff1_exprs)
#' ff2 <- flowFrame(ff2_exprs)
#' fs <- flowSet(ff1, ff2)
#' dim(ff1)
#'#>     events parameters
#'#>       1000         10
#' dim(ff2)
#'#>     events parameters
#'#>       1000         10
#' fs <- flowSet(ff1, ff2)
#' pooled_ff <- pool_fs(fs)
#' dim(pooled_ff)
#'#>     events parameters
#'#>       2000         10
#'
#' @export
pool_fs <- function(fs, copy_desc = TRUE) {
  # Check fs is a flowset ---------------------------
  stopifnot(grepl(pattern = "flowset", x = class(fs), ignore.case = TRUE))

  # Add cell ids as rownames ---------------------------
  add_cellids <- function(flow_frame, ff_label = NA) {
    if (is.na(ff_label) & exists("fs")) ff_label <- rownames(fs@phenoData)
    ff_dims <- dim(flowCore::exprs(flow_frame))
    ff_padded <- formatC(x = 1:ff_dims[[1]],
                         width = nchar(ff_dims[[1]]),
                         format = "d",
                         flag = "0")
    ff_ids <- base::paste(ff_label, ff_padded, sep = "_")
    rownames(flowCore::exprs(flow_frame)) <- ff_ids
    return(flow_frame)
  }
  fs <- foreach(i = seq_along(fs)) %do% add_cellids(fs[[i]])

  # Extract keywords from description slot ---------------------------
  extract_common_keywords <- function(fs) {
    fs_list <- foreach(i = seq_along(fs)) %do% fs[[i]]@description
    common_keywords <- base::Reduce(f = intersect,
                                    x = foreach(i = seq_along(fs_list)) %do% names(fs_list[[i]]))
    common_list <- foreach(i = seq_along(fs_list)) %do% fs_list[[i]][common_keywords]
    list_unnest <- function(list) {
      unnested_list <- foreach(i = seq_along(list)) %do% paste0(list[[i]], collapse = ", ")
      names(unnested_list) <- names(list)
      return(unnested_list)
    }
    common_unnested_list <- foreach(i = seq_along(common_list)) %do% list_unnest(common_list[[i]])
    common_names <- unique(foreach(i = seq_along(common_unnested_list), .combine = c) %do% names(common_unnested_list[[i]]))

    keyword_df <- as.data.frame(foreach(i = seq_along(common_unnested_list), .combine = rbind) %do% common_unnested_list[[i]][common_names])
    tmp <- foreach(i = seq_along(keyword_df)) %do% unique(keyword_df[, i])
    names(tmp) <- colnames(keyword_df)
    common_keywords <- names(tmp[foreach(i = seq_along(tmp), .combine = c) %do% length(tmp[[i]]) == 1])
    return(common_keywords)
  }

  # Pool into flowFrame ---------------------------
  if (copy_desc == TRUE) {
    fs_common_keywords <- extract_common_keywords(fs)
    desc <- fs[[1]]@description[fs_common_keywords]
    fs_pooled <- flowFrame(exprs = foreach(i = seq_along(fs), .combine = rbind) %do% fs[[i]]@exprs,
                           parameters = parameters(fs[[1]]),
                           description = desc)
    fs_pooled@description$`$TOT` <- dim(fs_pooled)[[1]]
  } else {
    fs_pooled <- flowFrame(exprs = foreach(i = seq_along(fs), .combine = rbind) %do% fs[[i]]@exprs,
                           parameters = parameters(fs[[1]]))
  }
  return(fs_pooled)
}

