#' Calculate ellipse for boral ordination
#'
#' @param metadf
#' @param ord
#' @param group
#' @param reflev
#' @param conf
#' @param cols
#'
#' @return
#' @export
#'
#' @examples
#' ## NOTE: As per the boral help, the values below MUST NOT be used in a real application;
#' ## they are only used here to make the examples run quick!!!
#' example_mcmc_control <- list(n.burnin = 10, n.iteration = 100,
#'                              n.thin = 1)
#' data("mite")
#' miteBoral <- boral(mite, family = "negative.binomial",
#' mcmc.control = example_mcmc_control,
#' lv.control = list(num.lv = 2),
#' row.eff = "fixed",
#' save.model = FALSE, calc.ics = FALSE)
#' par(mfrow = c(2, 2))


# Need to patch this up.

calcEllipse <- function(df, ord, group, reflev, conf = 0.95){

  scoreslvs <- data.frame(lvsplot(ord,
                       return.vals = TRUE,
                       biplot = FALSE)$scaled.lvs)
  names(scoreslvs) <- c("lvs1", "lvs2")
  scoresdf <- data.frame(df, scoreslvs)

  cols <-  c("lvs1", "lvs2")

  prehumanDF <- scoresdf[scoresdf[group] == reflev, ]

  mat <- cov.wt(prehumanDF[cols])

  tp <- sqrt(qchisq(conf, 2))

  xy <- veganCovEllipse(cov = mat$cov,
                        center = mat$center, tp)



  ret <- data.frame(xy, ellipseCenter1 = mat$center[[cols[1]]],
                    ellipseCenter2 = mat$center[[cols[2]]])
  return(ret)

}