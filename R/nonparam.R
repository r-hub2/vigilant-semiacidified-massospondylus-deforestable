
####
# Old R version
T_sq_nonpar_precomp <- function(means_1, means_2, nrow1, nrow2, cov1, cov2) {

  diff <- as.vector(means_1 - means_2)
  S <- as.matrix((1/(nrow1 + nrow2 - 2)) * (((nrow1 - 1) * cov1) + ((nrow2 - 1) * cov2)))
  res <- t(diff) %*% chol2inv(chol(S)) %*% diff
  res <- nrow1*nrow2/(nrow1+nrow2) * res
  as.numeric(res)
}


#### New version of Non-Param Model ####
# Training of the non-parametric model
#
# @examples
# library(deforestable)
#
# forestdir <- system.file('extdata/Forest/', package = "deforestable")
# Nonforestdir <- system.file('extdata/Non-forest/', package = "deforestable")
#
# NonParModel <- NonParamTrain(forestdir = forestdir,
#                              Nonforestdir = Nonforestdir)
# summary(NonParModel)
#
# # Read the target image
# tg_dir <- system.file('extdata/', package = "deforestable")
# test_image <- read_data_raster('smpl_1.jpeg', dir = tg_dir)
#
# res <- Nonparam_classifier(rastData=test_image, n_pts=10, Model=NonParModel, parallel=FALSE)
# res <- classify(data=test_image, Model=NonParModel, n_pts=10, parallel=FALSE, progress = 'text')
# jpeg::writeJPEG(image=res, target='NonPartest_im.jpeg')
#
# @export
NonParamTrain <- function(forestdir, Nonforestdir, Forest_list=NULL, Non_Forest_list=NULL){

  cl <- match.call()

  if (is.null(Forest_list))  Forest_list <- list.files(path=forestdir)
  if (is.null(Non_Forest_list)) Non_Forest_list <- list.files(path=Nonforestdir)

  # Combine image names into matrix with Forest and Not Forest labels
  dataset <- rbind(cbind(Im = Forest_list, label = rep("Forest", length(Forest_list))),
                   cbind(Im = Non_Forest_list, label = rep("Not Forest", length(Non_Forest_list))))


  #### Reading the data into lists ###############################################

  if(sum(dataset[,2]=="Forest")>1){
    f_read_dataset <- plyr::llply(as.list(dataset[dataset[,2]=="Forest",1]),
                                  read_data, dir = forestdir)
    rhfnames <- as.vector(dataset[dataset[,2]=="Forest",1])
  } else if (sum(dataset[,2]=="Forest")==1){
    f_read_dataset <- list(read_data(dataset[dataset[,2]=="Forest",][1], dir = forestdir))
    rhfnames <- vector(dataset[dataset[,2]=="Forest",][1])
  } else {
    f_read_dataset <- list()
    rhfnames <- vector()
  }

  if(sum(dataset[,2]=="Not Forest")>1){
    nf_read_dataset <- plyr::llply(as.list(dataset[dataset[,2]=="Not Forest",1]),
                                   read_data, dir = Nonforestdir)
    rhnfnames <- as.vector(dataset[dataset[,2]=="Not Forest",1])
  } else if (sum(dataset[,2]=="Not Forest")==1){
    nf_read_dataset <- list(read_data(dataset[dataset[,2]=="Not Forest",][1], dir = Nonforestdir))
    rhnfnames <- vector(dataset[dataset[,2]=="Not Forest",][1])
  } else {
    nf_read_dataset <- list()
    rhnfnames <- vector()
  }

  ################################################################################

  #### Forest processing ####

  # sample lengths
  v_f_ns <- vector(mode='numeric')
  v_f_ns <- c(v_f_ns, nrow(f_read_dataset[[1]]))

  # sample means
  f_means <- matrix(data=colMeans(f_read_dataset[[1]]),
                      nrow=1, ncol=3)

  # sample covariances
  f_cov_list <- list(cov(f_read_dataset[[1]]))

  for (i_f in 2:length(f_read_dataset)) {
    v_f_ns <- c(v_f_ns, nrow(f_read_dataset[[i_f]]))
    f_means <- rbind(f_means, colMeans(f_read_dataset[[i_f]]))
    f_cov_list <- c(
                    f_cov_list,
                    list(cov(f_read_dataset[[i_f]]))
                    )
  }


  #### Non-forest processing ####

  # sample lengths
  v_nf_ns <- vector(mode='numeric')
  v_nf_ns <- c(v_nf_ns, nrow(nf_read_dataset[[1]]))

  # sample means
  nf_means <- matrix(data=colMeans(nf_read_dataset[[1]]),
                     nrow=1, ncol=3)

  # sample covariances
  nf_cov_list <- list(cov(nf_read_dataset[[1]]))

  for (i_f in 2:length(nf_read_dataset)) {
    v_nf_ns <- c(v_nf_ns, nrow(nf_read_dataset[[i_f]]))
    nf_means <- rbind(nf_means, colMeans(nf_read_dataset[[i_f]]))
    nf_cov_list <- c(
      nf_cov_list,
      list(cov(nf_read_dataset[[i_f]]))
    )
  }


  ################################################################################


  #### Forest distances ####

  ds_f <- vector(mode='numeric')

  for (i_f in 1:length(f_read_dataset)) {

    ds_ind <- vector(mode='numeric')
    for (i_F in 1:length(f_read_dataset)) {
      res <- T_sq_nonpar_precomp_cpp(means_1 = f_means[i_F,], means_2 = f_means[i_f,],
                                     nrow1 = v_f_ns[i_F], nrow2 = v_f_ns[i_f],
                                     cov1 = f_cov_list[[i_F]], cov2 = f_cov_list[[i_f]])
      ds_ind <- c(ds_ind, res)
    }
    ds_f <- c(ds_f, min(ds_ind[ds_ind!=0]))
  }


  #### Non-forest distances ####

  ds_nf <- vector(mode='numeric')

  for (i_nf in 1:length(nf_read_dataset)) {

    ds_ind <- vector(mode='numeric')
    for (i_F in 1:length(f_read_dataset)) {
      res <- T_sq_nonpar_precomp_cpp(means_1 = f_means[i_F,], means_2 = nf_means[i_nf,],
                                     nrow1 = v_f_ns[i_F], nrow2 = v_nf_ns[i_nf],
                                     cov1 = f_cov_list[[i_F]], cov2 = nf_cov_list[[i_nf]])
      ds_ind <- c(ds_ind, res)
    }
    ds_nf <- c(ds_nf, min(ds_ind))
  }

  ll <- CompAccNP(ds_f, ds_nf)

  structure(list(call=cl,
                 acc=ll[[2]],
                 tp = ll[[3]],
                 fp = ll[[4]],
                 tn = ll[[5]],
                 fn = ll[[6]],
                 thres=ll[[1]],
                 forest_ns=v_f_ns,
                 forest_means=f_means,
                 forest_covars=f_cov_list
                ),
            class = c("ForestTrainNonParam", "ForestTrain"))
}



#### New classifier ####
# @export
Nonparam_classifier <- function(rastData, n_pts, Model, parallel=FALSE, progress = c("none","text")){

  thres <- as.numeric(Model$thres)
  matdata <- list(
                  as.matrix(rastData[[1]], wide=T),
                  as.matrix(rastData[[2]], wide=T),
                  as.matrix(rastData[[3]], wide=T)
                 )
  Gst <- Group_stats(matdata=matdata, n_pts=n_pts, fun=MultipleTester2, progress = progress,
                     parallel = parallel, paropts=list(.packages = c("deforestable")), Model=Model)
  Gst_clfd <- classifier_green(data=Gst, thres=thres)

  cc <- max(Gst$cc)
  data_rest <- clsfd_matr_restor(c=cc, n_pts=n_pts, data=Gst_clfd)
  data_rest
}



