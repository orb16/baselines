
#' Calculates distance from ellipse and distance from centroid.
#'
#' @param metadf metadata (e.g. groups). Needs to be same length and order as rows that went into ordination.
#' @param ord the ordination object. Currently only works with NMDS
#' @param group the column in the metadata which contains the grouping variable for the baseline
#' @param reflev the value in the column group which is the baseline level
#' @param ordiType what sort of ordiellipse? Options are all those which ordiellipse accepts
#' @param addConf use 95 CI or not? See ordiellipse options.
#'
#' @import sf
#' @import vegan
#' @import boral
#' @importFrom stats cov.wt qchisq
#' @return A named list: "distSF" with the distances from baseline as a spatial point object, "distDF", a data.frame version of the same, and "baseline_polygon" an SF object of the reference baseline polygon
#' @export
#'
#' @examples
#'
#' data("mite")
#' data("mite.env")
#' set.seed(999)
#' met <- metaMDS(mite, "jaccard", try = 50)
#' dlist <- calcEllipseDists(metadf = mite.env, ord = met,
#' group = "Topo", reflev = "Hummock")
#' par(mfrow = c(2, 2))
#' plot(met, type = "n")
#' plot(dlist[["baseline_polygon"]], add = TRUE,
#' col = adjustcolor("forestgreen", 0.2),
#' border = NA)
#' with(dlist[["distDF"]][dlist[["distDF"]]$Topo == "Hummock", ],
#' points(x = NMDS1, y = NMDS2, col = "forestgreen", pch = 16))
#' with(dlist[["distDF"]][dlist[["distDF"]]$Topo == "Blanket", ],
#' points(x = NMDS1, y = NMDS2, col = "black", pch = 16))
#' legend("topleft", pch = c(16, 16, 15),
#' col = c("forestgreen", "black",
#' adjustcolor("forestgreen", 0.5)),
#' legend = c("Hummock", "Blanket", "95% CI Ellipse"))
#' mtext(side = 3, "95% CI around centroid calculated")
#'
#'
#' plot(met, type = "n")
#' plot(dlist[["baseline_polygon"]], add = TRUE,
#' col = adjustcolor("forestgreen", 0.2),
#' border = NA)
#' mtext(side = 3, "Points sized by distance from baseline")
#' with(dlist[["distDF"]][
#' dlist[["distDF"]]$Topo == "Hummock", ],
#' points(x = NMDS1, y = NMDS2, cex = distEllipse + 0.5,
#' col = "forestgreen"))
#' with(dlist[["distDF"]][
#' dlist[["distDF"]]$Topo == "Blanket", ],
#' points(x = NMDS1, y = NMDS2, cex = distEllipse + 0.5,
#' col = "black"))
#'
#'
#' plot(met, type = "n")
#' plot(dlist[["baseline_polygon"]], add = TRUE,
#' col = adjustcolor("forestgreen", 0.2),
#' border = NA)
#' mtext(side = 3, "Blanket points by shrub prevalence")
#' colDF <- data.frame(Shrub = c("None", "Few", "Many"),
#' cols = I(c("skyblue", "cornflowerblue", "darkblue")))
#' with(dlist[["distDF"]][dlist[["distDF"]]$Topo == "Hummock", ],
#' points(x = NMDS1, y = NMDS2,  col = "grey"))
#' with(dlist[["distDF"]][
#' dlist[["distDF"]]$Topo == "Blanket", ],
#' points(x = NMDS1, y = NMDS2, col =  colDF[
#' match(Shrub, colDF$Shrub), "cols"]))
#' legend("topleft", pch = 1, legend = colDF$Shrub, col = colDF$cols)
#'
#' with(dlist[["distDF"]][dlist[["distDF"]]$Topo == "Blanket", ],
#' plot(x = Shrub, y = distEllipse, xlab = "Shrub",
#' ylab = "Distance from baseline"))
#' mtext(side = 3, "Distance from baseline as a function\nof another col in dataset (shrubs)")
#'
#' # Example with Boral
#' ## NOTE: As per the boral help, the values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100,
#'                              n.thin = 1)
#' data("mite")
#' data("mite.env")
#' miteBoral <- boral::boral(mite, family = "negative.binomial",
#' mcmc.control = example_mcmc_control,
#' lv.control = list(num.lv = 2),
#' row.eff = "fixed",
#' save.model = FALSE, calc.ics = FALSE)
#'
#' dlist <- calcEllipseDists(metadf = mite.env, ord = miteBoral,
#' group = "Topo", reflev = "Hummock")
#' par(mfrow = c(1, 2))
#' boral::lvsplot(miteBoral, biplot = FALSE,
#' col = "transparent",
#' main = "Boral ordination")
#' plot(dlist[["baseline_polygon"]], add = TRUE,
#' col = adjustcolor("forestgreen", 0.2),
#' border = NA)
##' plot(st_geometry(dlist[["distSF"]][dlist[["distSF"]]$Topo == "Hummock", ],),
#' add = TRUE,
#' col = "forestgreen",
#' pch = 16
#' )
#' plot(st_geometry(dlist[["distSF"]][dlist[["distSF"]]$Topo == "Blanket", ],),
#'      add = TRUE,
#'      col = "black",
#'      pch = 16
#' )
#' legend("topleft", pch = c(16, 16, 15),
#' col = c("forestgreen", "black",
#' adjustcolor("forestgreen", 0.5)),
#' legend = c("Hummock", "Blanket", "95% CI Ellipse"))
#' mtext(side = 3, "95% CI around centroid calculated")
#'
#' plot(met, type = "n", main = "NMDS ordination")
#' plot(dlist[["baseline_polygon"]], add = TRUE,
#' col = adjustcolor("forestgreen", 0.2),
#' border = NA)
#' plot(st_geometry(dlist[["distSF"]][dlist[["distSF"]]$Topo == "Hummock", ],),
#' add = TRUE,
#' col = "forestgreen",
#' pch = 16
#' )
#' plot(st_geometry(dlist[["distSF"]][dlist[["distSF"]]$Topo == "Blanket", ],),
#'      add = TRUE,
#'      col = "black",
#'      pch = 16
#' )
#' legend("topleft", pch = c(16, 16, 15),
#' col = c("forestgreen", "black",
#' adjustcolor("forestgreen", 0.5)),
#' legend = c("Hummock", "Blanket", "95% CI Ellipse"))
#' mtext(side = 3, "95% CI around centroid calculated")



calcEllipseDists <- function(metadf,
                        ord,
                        group = "Time_period",
                        reflev = "Pre-human",
                        ordiType = "sd", addConf = TRUE){

  if(!is.logical(addConf)){
    stop("addConf can only be TRUE or FALSE. You can't specify a numeric value sorry")
  }


  if(!ordiType %in% c("sd", "se")){
    stop("only types se and sd accepted")
  }

  # check that the column that includes the baseline is in the dataset
  if(! group %in% names(metadf)){
    stop(
      paste0("The group variable", " '", group, "' ", "is not found in the metadf object, the column names of which are ", "'", paste(names(metadf), collapse = ", "), "'")
      )
  }

  # check that the reference level in the 'group' column is present
  if(! reflev %in% as.character(unique(metadf[[group]]))){
    stop(
      paste0("reflev is specified as '", reflev, "'"), "but is not present in the group column", group)

  }

  # check the reference level is > 1 obs
  testColN <- sum(metadf[[group]] == reflev)
  if(! testColN > 1){
    stop("Reference group needs more than one obs")
  }

  # one method for vegan ordinations,
  # another method for boral
  if("metaMDS" %in% class(ord)){

    # calculate the scores
    metaScores <- data.frame(metadf, vegan::scores(ord, "sites"))


    cols <- c("NMDS1", "NMDS2")

    if(addConf){
      # this could be amended to take any number as per the vegan option
      # currently set to 0.95 and a true/false option
      lt <- vegan::ordiellipse(ord, group = metaScores[[group]],
                               kind = ordiType, draw = "none", conf = 0.95, display = "sites")
    } else {
      lt <- vegan::ordiellipse(ord, group = metaScores[[group]],
                               kind = ordiType, draw = "none", display = "sites")
    }

    refEllipseDF <- cbind(as.data.frame(with(metaScores[metaScores$group==reflev,],
                                             veganCovEllipse(lt[[reflev]]$cov,
                                                             lt[[reflev]]$center,
                                                             lt[[reflev]]$scale))),
                          group = reflev,
                          centroidNMDS1 =  lt[[reflev]]$center[[1]],
                          centroidNMDS2 = lt[[reflev]]$center[[2]])

    centreAxis1 <- lt[[reflev]]$center[[1]]
    centreAxis2 <- lt[[reflev]]$center[[2]]

  } else if("boral" %in% class(ord)){
    # where ordination is from boral
    ordScores <- data.frame(lvsplotData(ord,
                                        return.vals = TRUE,
                                        biplot = FALSE)$scaled.lvs)
    names(ordScores) <- c("lvs1", "lvs2")
    metaScores <- data.frame(metadf, ordScores)

    cols <-  c("lvs1", "lvs2")

    # prehuman ellipse
    prehumanDF <- metaScores[metaScores[group] == reflev, ]

    w <- rep(1, nrow(prehumanDF)) # we don't offer the user to set weights

    mat <- cov.wt(prehumanDF[cols], w)

    if(ordiType=="se"){
      mat$cov <- mat$cov * sum(mat$wt^2)
    }
    if (mat$n.obs == 1){
      mat$cov[] <- 0
    }


    if(addConf){
      tp <- sqrt(qchisq(0.95, 2))
    } else{
      tp <- 1    # this replicates the vegan package treatment
      # of the scale parameter
    }

    if (mat$n.obs > 1) {
      xy <- veganCovEllipse(cov = mat$cov, center = mat$center, tp)

    } else {
      stop("Need more than one point in reference group") # should be superfluous - tested at start
    }


    refEllipseDF <- data.frame(xy,
                               group = reflev,
                               ellipseCenterLVS1 = mat$center[[cols[1]]],
                               ellipseCenterLVS2 = mat$center[[cols[2]]])

    centreAxis1 <- mat$center[[cols[1]]]
    centreAxis2 <- mat$center[[cols[2]]]

  } else {
    stop("Only NMDS in vegan and boral ordinations supported. The ordination object is not either of these")
  }






  # calculate spatial points
  # NA crs operates in cartesian space, which is what we want
  spts <- sf::st_as_sf(metaScores,
                       coords = cols,
                       crs = NA)



  # reference polygon as sf object
  # takes and repeats the firsts point, to close the polygon
  # this is just the polygon
  phPoly <- st_polygon(list(as.matrix(refEllipseDF[c(1:nrow(refEllipseDF), 1), cols])))

  # proper sf object of reference level polygon
  phPolyData <- st_sfc(phPoly) %>%
    st_as_sf() %>%
    rename_geometry("geometry") %>%
    dplyr::mutate({{group}} := paste(reflev))

  # proper sf object of centroid point
  cent <- sf::st_as_sf(data.frame(dim1 =
                                    centreAxis1,
                                  dim2 = centreAxis2,
                                  row.names = NULL),
                       coords = c("dim1", "dim2"),
                       crs = NA)



  # calculate both sets of distances at once -
  # using the point data as the comparison
spts$distEllipse <- sf::st_distance(spts, phPolyData)
spts$disCentroid <- sf::st_distance(spts, cent)


dists <- spts %>% dplyr::mutate(centroid1 = .env$centreAxis1,
           centroid2 = .env$centreAxis2) %>%
    dplyr::select(names(spts %>% st_drop_geometry()), .data$centroid1, .data$centroid2,
                  .data$distEllipse, .data$disCentroid)

  if(ncol(st_distance(spts, phPolyData)) > 1){
    stop("error too many polygons")
  }


  if(ncol(st_distance(spts, cent)) > 1){
    stop("error too many centroids")
  }





  # to get the old data.frame form
  dfOnlyCoords <- st_coordinates(dists) %>% data.frame()
  names(dfOnlyCoords) <- cols
  dfOnly <- data.frame(dists %>% st_drop_geometry(), dfOnlyCoords)

  return(list(
    "distSF" = dists,
    "distDF" = dfOnly,
    "baseline_polygon" = phPolyData # this should be the sf polygon, but rename the polygon
  ))
}


