#' Caculate ecological distance from first sample
#'
#' @param speciesData dataset with species abundance, pre-transformed if desired
#' @param metaData  dataset with metadata, such as depth and age
#' @param idCol name of column in metaData by which speciesData is ordered. Generally depth or year. Needs to be in ascending order (ie getting deeper or older), in the same order as the species data.
#' @param idColType need to specify whether the idcol is a depth style or year style variable.
#' @param distMethod the distance method to use. See ?vegan::vegdist for options.
#' @param threshold keeps species only where those with a total summed abundance (not frequency) of more than or equal to the threshold. Default is zero.
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
#' idCol = "Age", distMethod = "jaccard", idColType = "year")
#'
#' with(distFirst, plot(x = Age, y = dist_jaccard,
#' ylab = "Jaccard distance"))
#'
calculateDistanceStart <- function(speciesData, metaData, idCol, idColType = NULL, distMethod,
                                   threshold = 0){

  # check that idCol is in the metadata df
  if(! idCol %in% names(metaData)){
    stop(
      cat("the idCol column must be found in the metaData dataset. idCol is",
             paste0("'", idCol, "'"), "the names in the dataset are",
          paste0("'", names(metaData), "'"))
      )
  }

  # check that the metaData is in order
  if(tolower(idColType) == "year"){
    if(! all.equal(rank(metaData[[idCol]]),  seq(from = 1, to = length(metaData[[idCol]]),
                                                 by = 1))) {
      stop("dataframe rows are out of order;\ndata must currently be formatted shallow to deep\nand in order. If reordering,\nmake sure you reorder the species data too")
    }
  } else if(tolower(idColType) == "depth"){
    if(! all.equal(rank(metaData[[idCol]]),  seq(to = 1, from = length(metaData[[idCol]]),
                                                 by = -1))) {
      stop("dataframe rows are out of order;\ndata must currently be formatted shallow to deep\nand in order. If reordering,\nmake sure you reorder the species data too")
    }
  } else {
    stop(paste("idColType must be either 'year' or 'depth'"))
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


  speciesData <- as.data.frame(speciesData[colSums(speciesData) >= threshold])


  metaData <- as.data.frame(metaData)

  row.names(speciesData) <- metaData[[idCol]]

  df <- as.data.frame(as.matrix(vegan::vegdist(speciesData, method = distMethod)))

  startPos <- which.min(names(df))

  startDist <- data.frame(df[1]); names(startDist) <- paste("dist", distMethod, sep = "_")
  startDist$id <- row.names(startDist)
  names(startDist)[grep("id", names(startDist))] <- idCol
  startDist[[idCol]] <- as.numeric(startDist[[idCol]])

  if(! all.equal(startDist[[idCol]], metaData[[idCol]])) {
    stop("Problem with non-equal id columns")
  }

  ret2 <- merge(startDist, metaData, by = idCol)
  ret <- ret2[c(names(metaData), paste("dist", distMethod, sep = "_"))]

  return(ret)

}

