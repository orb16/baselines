#' Calculate ecological distance from first sample
#'
#' Here we calculate ecological distance from the first sample in a dataset. Note that if you are calculating distance from the earliest sample in your dataset, this means your data will need to be reversed compared to the way it is typically ordered! i.e. if ordering by depth, you will ensure the 'deepest' sample is the first row in the dataframe. The other important note is that when you reorder your data, ensure that you reorder both the species and metadata (i.e. don't just reorder your metadata!).
#'
#' @param speciesData dataset with species abundance, pre-transformed if desired
#' @param metaData  dataset with metadata, such as depth and age
#' @param idCol name of column in metaData by which speciesData is ordered. Needs to be ordered such that the "start" sample is in the first row. Needs to be in the same order as the species data.
#' @param distMethod the distance method to use. See ?vegan::vegdist for options.
#' @param threshold keeps species only where those with a total summed abundance (not frequency) of more than or equal to the threshold. Default is zero. Note that if you transformed your data, the threshold is applied to the transformed data. This will likely affect what a reasonable threshold is.
#'
#' @return A dataframe with idCol and the distance for each point, including the distance from the first sample to itself. The function prints the value of the 'start' or reference level sample to which everything else is compared, and the name of the column it is derived from. This can be accessed using attributes(object)$refLevel if desired, where object is the name of the returned object.
#' @export
#'
#' @examples
#'
#' ## This example uses data from the analogue package
#' # note you should consider transforming your data if you haven't already
#' # the transformation step is not shown here
#'
#' require(analogue)
#' data(abernethy)
#'
#'# put the deepest sample first
#' abernethy2 <- dplyr::arrange(abernethy, -Depth)
#'
#' speciesData <- abernethy2[!names(abernethy2) %in% c("Depth", "Age")]
#' metaData <- abernethy2[c("Depth", "Age")]
#'
#' distFirst <- calculateDistanceStart(speciesData, metaData,
#' idCol = "Depth",
#'  distMethod = "jaccard")
#'
#' with(distFirst, plot(x = Age, y = dist_jaccard,
#' ylab = "Jaccard distance"))
#'
#' # access the reference level of the sample
#' attr(distFirst, "refLevel")
#'
calculateDistanceStart <- function(speciesData, metaData,
                                   idCol,
                                   distMethod,
                                   threshold = 0){



  # check that idCol is in the metadata df
  if(! idCol %in% names(metaData)){
    stop(
      cat("the idCol column must be found in the metaData dataset. idCol is",
          paste0("'", idCol, "'"), "the names in the dataset are",
          paste0("'", names(metaData), "'"))
    )
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


  # warn on threshold
  numCol <- ncol(speciesData)
  numSpec <- ncol(speciesData[colSums(speciesData) > 0])
  numThresh <- ncol(speciesData[colSums(speciesData) >= threshold])
  if(numCol > numSpec){
    warning("Some summed species abundances total zero")
  }

  if(numThresh < 1){
    stop("No species meet the specified threshold! Please reconsider the threshold value")
  }

  if(numThresh < numSpec){
    warning(paste("Original data contained", numSpec, "species with a summed abundance of >0. After applying the threshold,", numThresh, "species remain"))
  }

  # filter such that species included only if above the threshold (default = 0)
  speciesData <- as.data.frame(speciesData[colSums(speciesData) >= threshold])

  metaData <- as.data.frame(metaData)

  # apply rownames to species data
  # this is why the species and meta data have to be in the same order
  row.names(speciesData) <- metaData[[idCol]]

  # calculate distance, convert to dataframe
  df <- as.data.frame(as.matrix(vegan::vegdist(speciesData, method = distMethod)))


  # select the first column (which is why it is important to have) the first row of the df being the "start" sample
  startDist <- data.frame(df[1])

  # record the name of the distance on the distance var
  names(startDist) <- paste("dist", distMethod, sep = "_")

  #convert the names from the distance matrix to an explicit column
  startDist$id <- row.names(startDist)

  # then rename to the original name
  names(startDist)[grepl("^id$", names(startDist))] <- idCol

  # as.character conversion here conservative - unlikely to be required
  startDist[[idCol]] <- as.numeric(as.character(startDist[[idCol]]))

  # a hopefully redundant test
  if(! identical(startDist[[idCol]], metaData[[idCol]])) {
    stop("Problem with non-equal id columns")
  }

  # join the distance data back to the metadata
  ret <- merge(metaData, startDist, by = idCol,
               sort = FALSE)



  # record the reference level as an attribute
  attr(ret, "refLevel") <- paste0("Start/reference level value is ", names(df[1]),
                                     " from column '", idCol, "'")

  # print the reference level to remind folks
  print(attr(ret, "refLevel"))
  return(ret)


}

