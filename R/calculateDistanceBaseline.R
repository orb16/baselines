
#' Calculates distance from ellipse and distance from centroid.
#'
#' @param metadf metadata (e.g. groups). Needs to be same length and order as rows that went into ordination.
#' @param ord the ordination object. Currently only works with NMDS
#' @param group the column in the metadata which contains the grouping variable for the baseline
#' @param reflev the value in the column group which is the baseline level
#' @param ordiType what sort of ordiellipse? Options are all those which ordiellipse accepts
#' @param conf use 95 CI or not? See ordiellipse options.
#'
#' @import sp
#' @import rgeos
#' @import vegan
#' @return A list, with a dataframe containing the metadf and distance from baseline, second, a spatial points object. Third, the baseline ellipse created as a spatial object.
#' @export
#'
#' @examples
#'
#' data("mite")
#' data("mite.env")
#' met <- metaMDS(mite, "jaccard")
#' dlist <- calcMyDists(metadf = mite.env, ord = met, group = "Topo", reflev = "Hummock")
#' par(mfrow = c(1, 3))
#' plot(met, type = "n")
#' plot(dlist[["baseline_polygon"]], add = TRUE,
#' col = adjustcolor("forestgreen", 0.2),
#' border = NA)
#' points(dlist[["all_points"]][dlist[["all_points"]]$Topo == "Hummock", ], col = "forestgreen",
#' pch = 16)
#' points(dlist[["all_points"]][dlist[["all_points"]]$Topo == "Blanket", ], col = "black",
#' pch = 16)
#' with(dlist[["distDF"]][dlist[["distDF"]]$Topo == "Blanket", ],
#' plot(x = Shrub, y = distEllipse, xlab = "Shrub",
#' ylab = "Distance from baseline"))
#' plot(met, type = "n")
#' plot(dlist[["baseline_polygon"]], add = TRUE,
#' col = adjustcolor("forestgreen", 0.2),
#' border = NA)
#' with(dlist[["distDF"]][dlist[["distDF"]]$Topo == "Hummock", ],
#' points(x = NMDS1, y = NMDS2, cex = distEllipse, col = "forestgreen"))
#' with(dlist[["distDF"]][dlist[["distDF"]]$Topo == "Blanket", ],
#' points(x = NMDS1, y = NMDS2, cex = distEllipse, col = "black"))

calcMyDists <- function(metadf,
                        ord,
                        group = "Time_period",
                        reflev = "Pre-human",
                        ordiType = "sd", conf = TRUE){

  # if (!requireNamespace("sp", quietly = TRUE)) {
  #   stop("Package \"sp\" needed for this function to work. Please install it.",
  #        call. = FALSE)
  # }
  # if (!requireNamespace("rgeos", quietly = TRUE)) {
  #   stop("Package \"rgeos\" needed for this function to work. Please install it.",
  #        call. = FALSE)
  # }

  metaScores <- data.frame(metadf, vegan::scores(ord))

  allPts2 <- metaScores[, c("NMDS1", "NMDS2")]
  spts <- sp::SpatialPoints(allPts2)
  sptsDF <- sp::SpatialPointsDataFrame(spts,
                                       data = metaScores[! grepl("NMDS1|NMDS2", names(metaScores))])

  if(conf){
    lt <- vegan::ordiellipse(ord, group = metaScores[[group]],
                      kind = ordiType, draw = "none")
  } else {
    lt <- vegan::ordiellipse(ord, group = metaScores[[group]],
                      kind = ordiType, draw = "none", conf = 0.95)
  }


  def_ell <- data.frame() #sets up a data frame before running the function.
  for(g in c(reflev)){

    tmp <- cbind(as.data.frame(with(metaScores[metaScores$group==g,],
                                    veganCovEllipse(lt[[g]]$cov,lt[[g]]$center,lt[[g]]$scale)))
                 ,group=g)
    def_ell <- rbind(def_ell, tmp)

    if(g == reflev){


      phPoly <- sp::SpatialPolygons(list(sp::Polygons(srl =
                                                        list(sp::Polygon(as.matrix(tmp[1:2]), hole = FALSE)), ID = 1)))

      phPolyData <- sp::SpatialPolygonsDataFrame(phPoly,
                                                 data = data.frame(lev = reflev))
      names(phPolyData@data) <- paste(group)



      dists <- apply(rgeos::gDistance(spts, phPoly, byid=TRUE),2,min)

      cent <- sp::SpatialPoints(data.frame(NMDS1 = lt[[g]]$center[1], NMDS2 = lt[[g]]$center[2], row.names = NULL))

      centDist <- rgeos::gDistance(spts, cent, byid = TRUE)

      distPH <- data.frame(metaScores, distEllipse = dists,
                           disCentroid = as.vector(centDist),
                           centroidNMDS1 =  lt[[g]]$center[1],
                           centroidNMDS2 = lt[[g]]$center[2])

    }

  }
    return(list("distDF" = distPH, "baseline_polygon" = phPolyData, "all_points" = sptsDF))
  }
