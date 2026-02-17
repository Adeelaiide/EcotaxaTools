#' dendrocut_auto
#'
#'Estimate the best cut in a hierarchical cluster analysis based on the 
#'inverted cumulative sum of the coefficient of inconsistance.
#'
#'@param hc a cluster list obtain from hclust
#'
#' @return A numerical value of the cut
#' @export
#'
#' @examples dendrocut_auto(hclust(vegdist(data_stat, method="euclidean"), method="complete"))
#' 
dendrocut_auto <- function(hc){
  
  heights <- hc$height
  
  #Inconsistency coefficient approximation
  incon <- matrix(NA, nrow = length(heights), ncol = 4)
  for(i in 1:length(heights)){
    lower_heights <- heights[1:i]
    incon[i,1] <- mean(lower_heights)
    incon[i,2] <- ifelse(i==1, 0, sd(lower_heights))
    incon[i,3] <- length(lower_heights)
    incon[i,4] <- ifelse(incon[i,2] == 0, 0,
                         (heights[i] - incon[i,1]) / incon[i,2])
  }
  
  # Inverted cumulative sum
  test <- cumsum(rev(incon[,4]))
  I <- which(test == max(test)) #Maximum index
  dendrocut_auto <- incon[I[1],1] #Suggested value
  return(dendrocut_auto)
}