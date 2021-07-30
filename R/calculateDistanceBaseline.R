
#' Calculates distance from ellipse and distance from centroid.
#'
#' @param metadf metadata (e.g. groups). Needs to be same length and order as rows that went into ordination.
#' @param ord the ordination object. Currently only works with NMDS
#' @param group the column in the metadata which contains the grouping variable for the baseline
#' @param reflev the value in the column group which is the baseline level
#' @param ordiType what sort of ordiellipse? Options are all those which ordiellipse accepts
#' @param addConf use 95 CI or not? See ordiellipse options.
#'
#' @import sp
#' @import rgeos
#' @import vegan
#' @import boral
#' @importFrom stats cov.wt qchisq
#' @return A list, with a dataframe containing the metadf and distance from baseline, second, a spatial points object. Third, the baseline ellipse created as a spatial object.
#' @export
#'
#' @examples
#'
#' set.seed(999)
#' data("mite")
#' data("mite.env")
#' met <- metaMDS(mite, "jaccard", try = 100)
#' dlist <- calcEllipseDists(metadf = mite.env, ord = met,
#' group = "Topo", reflev = "Hummock")
#' par(mfrow = c(2, 2))
#' plot(met, type = "n")
#' plot(dlist[["baseline_polygon"]], add = TRUE,
#' col = adjustcolor("forestgreen", 0.2),
#' border = NA)
#' points(dlist[["all_points"]][
#' dlist[["all_points"]]$Topo == "Hummock", ],
#' col = "forestgreen",
#' pch = 16)
#' points(dlist[["all_points"]][
#' dlist[["all_points"]]$Topo == "Blanket", ],
#' col = "black",
#' pch = 16)
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
#' points(dlist[["all_points"]][dlist[["all_points"]]$Topo == "Hummock", ],
#' col = "forestgreen",
#' pch = 16)
#' points(dlist[["all_points"]][dlist[["all_points"]]$Topo == "Blanket", ],
#' col = "black",
#' pch = 16)
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
#' points(dlist[["all_points"]][dlist[["all_points"]]$Topo == "Hummock", ], col = "forestgreen",
#' pch = 16)
#' points(dlist[["all_points"]][dlist[["all_points"]]$Topo == "Blanket", ], col = "black",
#' pch = 16)
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

  if(! group %in% names(metadf)){
    stop(
      cat("The column name", paste0("'", group, "'"), "is not found in the metadf object, the column names of which are", paste0("'", names(metadf), "'"))
      )
  }

  if(! reflev %in% as.character(unique(metadf[[group]]))){
    stop(
      cat(paste0("reflev '", reflev, "'"), "needs to be the column", group)
         )
  }

  # one method for vegan ordinations,
  # another method for boral
  if("metaMDS" %in% class(ord)){

    # calculate the scores
    metaScores <- data.frame(metadf, vegan::scores(ord, "sites"))


    cols <- c("NMDS1", "NMDS2")

    if(addConf){
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

    ordScores <- data.frame(lvsplotData(ord,
                                    return.vals = TRUE,
                                    biplot = FALSE)$scaled.lvs)
    names(ordScores) <- c("lvs1", "lvs2")
    metaScores <- data.frame(metadf, ordScores)

    cols <-  c("lvs1", "lvs2")

    # prehauman ellipse
    prehumanDF <- metaScores[metaScores[group] == reflev, ]
    mat <- cov.wt(prehumanDF[cols])
    if(addConf){
      tp <- sqrt(qchisq(0.95, 2))
    } else{
      stop("Ellipse without 95% CI not yet implemented for boral")
    }

    xy <- veganCovEllipse(cov = mat$cov,
                          center = mat$center, tp)
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

  allPts2 <- metaScores[cols]
  # this can be done for both
  spts <- sp::SpatialPoints(allPts2)
  sptsDF <- sp::SpatialPointsDataFrame(spts,
                                       data = metaScores[! names(metaScores) %in% cols])






  phPoly <- sp::SpatialPolygons(list(sp::Polygons(srl =
                                                    list(sp::Polygon(as.matrix(refEllipseDF[cols]), hole = FALSE)), ID = 1)))

  phPolyData <- sp::SpatialPolygonsDataFrame(phPoly,
                                             data = data.frame(lev = reflev))
  names(phPolyData@data) <- paste(group)



  dists <- apply(rgeos::gDistance(spts, phPoly, byid=TRUE), 2, min)

  # fix this
  cent <- sp::SpatialPoints(data.frame(dim1 = centreAxis1,
                                       dim2 = centreAxis2,
                                       row.names = NULL))

  centDist <- rgeos::gDistance(spts, cent, byid = TRUE)

  distPH <- data.frame(metaScores,
                       centroid1 = centreAxis1,
                       centroid2 = centreAxis2,
                       distEllipse = dists,
                       disCentroid = as.vector(centDist))



    return(list("distDF" = distPH,
                "baseline_polygon" = phPolyData,
                "all_points" = sptsDF,
                "baseline_polygon_DF" = refEllipseDF))
  }

