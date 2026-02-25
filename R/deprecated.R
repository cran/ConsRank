#' Defunct functions in ConsRank
#' 
#' These functions have been removed from ConsRank. Use the alternatives listed below.
#' 
#' @details
#' \itemize{
#'   \item \code{EMCons}: Use \code{\link{consrank}} with \code{algorithm = "BB"}
#'   \item \code{QuickCons}: Use \code{\link{consrank}} with \code{algorithm = "quick"}
#'   \item \code{BBFULL}: Use \code{\link{consrank}} with \code{algorithm = "BB"} and \code{full = TRUE}
#'   \item \code{FASTcons}: Use \code{\link{consrank}} with \code{algorithm = "fast"}
#'   \item \code{DECOR}: Use \code{\link{consrank}} with \code{algorithm = "decor"}
#'   \item \code{FASTDECOR}: Use \code{\link{consrank}} with \code{algorithm = "fastdecor"}
#'   \item \code{labels}: Functionality merged into \code{\link{rank2order}}
#' }
#' 
#' @name ConsRank-defunct
#' @keywords internal
NULL

#' @rdname ConsRank-defunct
#' @export
EMCons <- function(...) {
  .Defunct(
    new = "consrank",
    package = "ConsRank",
    msg = "EMCons() is defunct. Use consrank(algorithm = 'EM') instead."
  )
}

#' @rdname ConsRank-defunct
#' @export
QuickCons <- function(...) {
  .Defunct(
    new = "consrank",
    package = "ConsRank",
    msg = "QuickCons() is defunct. Use consrank(algorithm = 'quick') instead."
  )
}

#' @rdname ConsRank-defunct
#' @export
BBFULL <- function(...) {
  .Defunct(
    new = "consrank",
    package = "ConsRank",
    msg = "BBFULL() is defunct. Use consrank(algorithm = 'BB', full = TRUE) instead."
  )
}

#' @rdname ConsRank-defunct
#' @export
FASTcons <- function(...) {
  .Defunct(
    new = "consrank",
    package = "ConsRank",
    msg = "FASTcons() is defunct. Use consrank(algorithm = 'fast') instead."
  )
}

#' @rdname ConsRank-defunct
#' @export
DECOR <- function(...) {
  .Defunct(
    new = "consrank",
    package = "ConsRank",
    msg = "DECOR() is defunct. Use consrank(algorithm = 'decor') instead."
  )
}

#' @rdname ConsRank-defunct
#' @export
FASTDECOR <- function(...) {
  .Defunct(
    new = "consrank",
    package = "ConsRank",
    msg = "FASTDECOR() is defunct. Use consrank(algorithm = 'fastdecor') instead."
  )
}

#' @rdname ConsRank-defunct
#' @export
labels <- function(...) {
  .Defunct(
    new = "consrank",
    package = "ConsRank",
    msg = "labels() is defunct. This functionality is now integrated in rank2order()."
  )
}