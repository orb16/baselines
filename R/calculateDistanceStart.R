#' Caculate ecological distance from first sample
#'
#' @param speciesData dataset with species abundance, pre-transformed if desired
#' @param metaData  dataset with metadata, such as depth and age
#' @param idCol name of column in metaData by which speciesData is ordered. Generally depth or year. Needs to be in ascending order (ie getting deeper or older), in the same order as the species data.
#' @param distMethod the distance method to use. See ?vegan::vegdist for options.
#'
#' @return A dataframe with idCol and the distance for each point, including the distance from the first sample to itself.
#' @export
#'
#' @examples
#'
#' ## This example uses data from the analogue package
#'
#' require(analogue)
#' data(abernethy)
#' speciesData <- abernethy[!names(abernethy) %in% c("Depth", "Age")]
#' metaData <- abernethy[c("Depth", "Age")]
#'
#' distFirst <- calculateDistanceStart(speciesData, metaData,
#' idCol = "Age", distMethod = "jaccard")
#'
#' with(distFirst, plot(x = Age, y = dist_jaccard,
#' ylab = "Jaccard distance"))
#'
calculateDistanceStart <- function(speciesData, metaData, idCol, distMethod){

  # check that idCol is in the metadata df
  if(! idCol %in% names(metaData)){
    stop(
      cat("the idCol column must be found in the metaData dataset. idCol is",
             paste0("'", idCol, "'"), "the names in the dataset are",
          paste0("'", names(metaData), "'"))
      )
  }

  # check that the metaData is in order
  if(! all.equal(rank(metaData[[idCol]]),  seq(from = 1, to = length(metaData[[idCol]]),
                                               by = 1))) {
    stop("dataframe rows are out of order;\ndata must currently be formatted shallow to deep\nand in order. If reordering,\nmake sure you reorder the species data too")
  }

  # make sure the distance method can be applied using vegdist
  if(!distMethod %in% c("manhattan", "euclidean", "canberra",
                        "clark", "bray", "kulczynski",
                        "jaccard", "gower", "altGower",
                        "morisita", "horn", "mountford",
                        "raup", "binomial", "chao",
                        "cao", "mahalanobis")){
    stop("Only distance measures available in ?vegan::vegdist can be used.")
  }



  row.names(speciesData) <- metaData[[idCol]]

  df <- as.data.frame(as.matrix(vegan::vegdist(speciesData, method = distMethod)))

  startPos <- which.min(names(df))

  startDist <- df[1]; names(startDist) <- paste("dist", distMethod, sep = "_")
  startDist$id <- row.names(startDist)
  names(startDist)[grep("id", names(startDist))] <- idCol

  ret <- startDist[c(2,1)]; row.names(ret) <- NULL

  return(ret)

}
