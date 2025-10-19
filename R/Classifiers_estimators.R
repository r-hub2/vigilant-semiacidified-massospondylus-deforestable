
#' Koutrouvelis parameter estimation of image data
#'
#' In data, there are three columns and each column corresponds to the color
#' intensity of one channel: red, green and blue correspondingly. The four
#' parameters: alpha, beta, gamma and delta, of the stable distribution is
#' estimated for each of these channels using the Koutrouvelis regressions-type
#' technique.
#' @param data matrix or data frame with color intensities of red, green and
#'   blue for an image.
#' @return a data frame with columns alpha, beta, gamma, delta and rows red,
#'   green and blue.
#' @examples
#' library(deforestable)
#'
#' Forestdir <- system.file('extdata/Forest/', package = "deforestable")
#' test_image <- read_data('_6_33_.jpeg', dir = Forestdir)
#'
#' pars <- Koutparams(test_image)
#'
#' pars
#'
#' @export
Koutparams <- function(data){

  r_m <- mean(data[,1])
  g_m <- mean(data[,2])
  b_m <- mean(data[,3])

  theta0r <- c(1.72, 0, 0.04, r_m)
  theta0g <- c(1.65, 0, 0.04, g_m)
  theta0b <- c(1.81, 0, 0.038, b_m)

  red <- KoutParametersEstim(x=data[,1], theta0 = NULL, spacing = "Kout",
                             pm = 1, tol = 0.05, NbIter = 10, PrintTime = FALSE)$Estim$par

  red[4] <- r_m

  green <-  KoutParametersEstim(x=data[,2], theta0 = NULL, spacing = "Kout",
                                pm = 1, tol = 0.05, NbIter = 10, PrintTime = FALSE)$Estim$par

  green[4] <- g_m

  blue <-  KoutParametersEstim(x=data[,3], theta0 = NULL, spacing = "Kout",
                               pm = 1, tol = 0.05, NbIter = 10, PrintTime = FALSE)$Estim$par

  blue[4] <- b_m

  dd <- rbind(red,green,blue)
  colnames(dd) <- c('alpha', 'beta', 'gamma', 'delta')
  as.data.frame(dd)
}



Mahala_dist <- function(params, dataset) {

  diff <- as.vector(colMeans(params) - colMeans(dataset[,1:3]))
  cov1 <- cov(params)
  cov2 <- cov(dataset[,1:3])
  nrow1 <- nrow(params)
  nrow2 <- nrow(dataset[,1:3])
  S <- as.matrix((1/(nrow1 + nrow2 - 2)) * (((nrow1 - 1) * cov1) + ((nrow2 - 1) * cov2)))
  res <- t(diff) %*% chol2inv(chol(S)) %*% diff
  res[[1]]
}



Group_stats <- function(matdata, n_pts, fun,
                        progress = "text", parallel = TRUE, paropts=NULL, ...){

  # number of groups
  r <- nrow(matdata[[1]]) %/% n_pts
  c <- ncol(matdata[[1]]) %/% n_pts

  # capture only one channel, small size
  m_df <- as.array(matdata)
  data_ch1 <- m_df[[1]]
  data_ch2 <- m_df[[2]]
  data_ch3 <- m_df[[3]]

  # Data truncation
  data_ch1 <- data_ch1[1:(r*n_pts), 1:(c*n_pts)]
  data_ch2 <- data_ch2[1:(r*n_pts), 1:(c*n_pts)]
  data_ch3 <- data_ch3[1:(r*n_pts), 1:(c*n_pts)]

  # Vectorization
  data_ch1_vec <- as.vector(data_ch1)
  data_ch2_vec <- as.vector(data_ch2)
  data_ch3_vec <- as.vector(data_ch3)

  # indexes of groups
  r1 <- rep(1:r, each = n_pts)
  rr <- rep(r1, times = c*n_pts)

  cc <- rep(1:c, each=r*n_pts^2)

  # data with indexes
  data_merged <- cbind(data_ch1_vec,
                       data_ch2_vec,
                       data_ch3_vec,
                       rr, cc)

  data_merged <- as.data.frame(data_merged)

  # result
  splitter_plyr <- utils::getFromNamespace("splitter_d", "plyr")
  data_split <- splitter_plyr(data_merged, plyr::.(rr, cc), drop = FALSE)
  data_l <- plyr::ldply(.data=data_split, .fun=fun, .paropts=paropts,
                        .progress = progress, .parallel = parallel, .inform=FALSE, ...)

  data_l

}


classifier_green <- function(data, thres){

  data.frame(
    cbind(data[,1:2], as.numeric(data[,3] < thres)),
    fix.empty.names=FALSE)
}


classifier_all <- function(data, thres){

  data.frame(
    cbind(data[,1:2], as.numeric(data[,3] < thres[1] & data[,4] < thres[2] & data[,5] < thres[3])),
    fix.empty.names=FALSE)
}


clsfd_matr_restor <- function(c, n_pts, data){

  rr <- NULL # avoids NOTEs when being built
  splitter_plyr <- utils::getFromNamespace("splitter_d", "plyr")

  test <- data[,1:3]
  test2 <- splitter_plyr(test, plyr::.(rr), drop = FALSE)


  test3 <- lapply(test2, function(x) x[,3])
  test3 <- lapply(test3, function(x) rep(x, times=n_pts, each=n_pts))

  test4 <- lapply(test3, function(x) matrix(x, nrow = n_pts, ncol = c*n_pts, byrow = TRUE))
  test4 <- do.call(rbind, test4)
  test4
}



# Non-parametric tester
MultipleTester2 <- function(data, Model){

  data <- data[,1:3]

  nrow_1 <- nrow(data)
  means_1 <- as.vector(colMeans(data))
  cov_1 <- cov(data)

  N_forests <- length(Model$forest_ns)

  Mah_dists <- vector(mode = 'numeric')
  for (ind in 1:N_forests) {
    Mah_dists <- c(Mah_dists,
                   T_sq_nonpar_precomp_cpp(means_1 = means_1, means_2 = as.vector(Model$forest_means[ind,]),
                                           nrow1 = nrow_1, nrow2 = Model$forest_ns[ind],
                                           cov1 = cov_1, cov2 = Model$forest_covars[[ind]]))
  }

  min(Mah_dists)
}


# Parametric tester
MultipleTesterParam2 <- function(data, Model){

  data <- data[,1:3]
  nrow_1 <- nrow(data)

  ECF_dists <- vector(mode = 'numeric')
  N_forests <- length(Model$f_Z)

  for (ind in 1:N_forests) {

    Z_1 <- Z_n(t_par=Model$t_pars[ind], X=data[,2])
    Z_2 <- Model$f_Z[[ind]]
    mtrx <- Model$inv_f_sigmas[[ind]]
    ECF_dists <- c(ECF_dists,
                   nrow_1 * as.numeric(t(Z_1-Z_2) %*% mtrx %*% (Z_1-Z_2)))
  }

  min(ECF_dists)
}


