#' Save a set of masks as a table
#'
#' Writes all masks of a set into one tab-separated file with three columns:
#' \emph{mask} (ordinal number), \emph{sample} and \emph{marker}. The sample
#' and marker names are taken from the dimnames of the masks when present,
#' i.e. when the masks were generated with the \emph{samples} and
#' \emph{markers} arguments of \code{\link{GenerateMaskSet}} or annotated
#' otherwise; ordinal numbers are written for the unannotated dimensions.
#' The table makes the hidden genotypes portable outside R.
#'
#' @param masks list of masks as binary matrices
#' @param file path/to/filename.tsv
#' @param samples optional character vector of length \emph{N} with sample
#'   ids, used for the rows of masks without row names
#' @param markers optional character vector of length \emph{M} with marker
#'   ids, used for the columns of masks without column names
#'
#' @return No values; the number of rows written is reported as a message
#' @export
#'
WriteMaskSet <- function(masks, file, samples = NULL, markers = NULL) {

  rows <- lapply(seq_along(masks), function(i) {

    m <- masks[[i]]
    if (!is.matrix(m) || any(dim(m) == 0)) return(NULL)

    s <- rownames(m)
    if (is.null(s)) s <- samples
    if (is.null(s)) s <- as.character(seq_len(nrow(m)))

    k <- colnames(m)
    if (is.null(k)) k <- markers
    if (is.null(k)) k <- as.character(seq_len(ncol(m)))

    idx <- which(m == 1, arr.ind = TRUE)
    if (nrow(idx) == 0) return(NULL)

    data.frame(mask = i,
               sample = s[idx[, "row"]],
               marker = k[idx[, "col"]],
               stringsAsFactors = FALSE)
  })

  tab <- do.call(rbind, rows)
  if (is.null(tab)) stop("The masks are empty, nothing to write", call. = FALSE)

  utils::write.table(tab, file, sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("%s hidden genotypes of %s masks are written to %s",
                  nrow(tab), length(masks), file))
}
