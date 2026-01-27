#' bv_to_esdum
#'
#' Compute ESD in micron from Biovolume in mm3
#'
#' @param x Biovolume in mm3
#'
#' @return ESD in micron
#' @export
#'
#' @examples bv_to_esdum(x=Biovolume)
bv_to_esdum <- function(x) {
  x<- (6 *x/pi)^(1/3)*1000
}

