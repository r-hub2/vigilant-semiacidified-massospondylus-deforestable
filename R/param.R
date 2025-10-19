
#### New version of Param Model ####
# Training of the non-parametric model
#
# @examples
# library(deforestable)
#
# forestdir <- system.file('extdata/Forest/', package = "deforestable")
# Nonforestdir <- system.file('extdata/Non-forest/', package = "deforestable")
#
# ParModel <- ParamTrain(forestdir = forestdir,
#                        Nonforestdir = Nonforestdir)
#
# # Read the target image
# tg_dir <- system.file('extdata/', package = "deforestable")
# test_image <- read_data_raster('smpl_1.jpeg', dir = tg_dir)
#
# res <- Param_classifier(rastData=test_image, n_pts=10, progress='text',
#                         Model=ParModel, parallel=FALSE)
# res <- classify(data=test_image, Model=ParModel,
#                 n_pts=10, parallel=FALSE, progress = 'text')
# jpeg::writeJPEG(image=res, target='Partest_im.jpeg')
# @export
ParamTrain <- function(forestdir, Nonforestdir, Forest_list=NULL, Non_Forest_list=NULL){

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
  # only the green channel is included

  v_f_ns <- vector(mode='numeric')
  t_pars <- vector(mode='numeric')
  z_n <- list()
  f_Z <- list()
  f_stab_params <- list()
  f_sigmas <- list()
  inv_f_sigmas <- list()

  for (ind in 1:length(f_read_dataset)) {

    # data lengths
    v_f_ns <- c(v_f_ns, nrow(f_read_dataset[[ind]]))

    # tables of stable parameters
    pars <- Koutparams(f_read_dataset[[ind]])
    f_stab_params <- c(f_stab_params, list(pars))
    t_par <- 1/15/pars['green','gamma']
    t_pars <- c(t_pars, t_par)

    # vectors of cfs
    z_n <- Z_n(t_par=t_par,
               X=f_read_dataset[[ind]][,'green']) 
    f_Z <- c(f_Z, list(z_n))

    # Cf covariance matrices
    sigma <- stbl_param_covmtrx_cpp(t_par=t_par, theta=base::as.vector(pars['green',], mode = 'numeric'))
    inv_f_sigma<- solve(sigma)
    f_sigmas <- c(f_sigmas, list(sigma))
    inv_f_sigmas <- c(inv_f_sigmas, list(inv_f_sigma))
  }

  #### Non-forest processing ####
  # only the green channel is included

  v_nf_ns <- vector(mode='numeric')
  # z_n <- list()
  nf_Z <- list()

  for (ind in 1:length(nf_read_dataset)) {
    # data lengths
    v_nf_ns <- c(v_nf_ns, nrow(nf_read_dataset[[ind]]))
  }


  ################################################################################


  #### Forest distances ####
  ds_f <- vector(mode='numeric')

  for (i_f in 1:length(f_read_dataset)) {
    ds_ind <- vector(mode='numeric')
    for (i_F in 1:length(f_read_dataset)) {

      z_n <- Z_n(t_par=t_pars[i_F],
                 X=f_read_dataset[[i_f]][,'green'])

      res <- v_f_ns[i_f] %*% t(z_n - f_Z[[i_F]]) %*% inv_f_sigmas[[i_F]] %*% (z_n - f_Z[[i_F]])
      res <- as.numeric(res)
      ds_ind <- c(ds_ind, res)
    }
    ds_f <- c(ds_f, min(ds_ind[ds_ind!=0]))
  }

  #### Non-forest distances ####

  ds_nf <- vector(mode='numeric')

  for (i_nf in 1:length(nf_read_dataset)) {
    ds_ind <- vector(mode='numeric')
    for (i_F in 1:length(f_read_dataset)) {

      z_n <- Z_n(t_par=t_pars[i_F],
                 X=nf_read_dataset[[i_nf]][,'green'])

      res <- v_nf_ns[i_nf] %*% t(z_n - f_Z[[i_F]]) %*% inv_f_sigmas[[i_F]] %*% (z_n - f_Z[[i_F]])
      res <- as.numeric(res)
      ds_ind <- c(ds_ind, res)
    }
    ds_nf <- c(ds_nf, min(ds_ind))
  }


  ####  ####

  ll <- CompAccNP(ds_f, ds_nf)

  structure(list(call=cl,
                 acc=ll[[2]],
                 tp = ll[[3]],
                 fp = ll[[4]],
                 tn = ll[[5]],
                 fn = ll[[6]],
                 thres=ll[[1]],
                 f_stab_params=f_stab_params,
                 forest_sigmas=f_sigmas,
                 inv_f_sigmas=inv_f_sigmas,
                 f_Z=f_Z,
                 t_pars=t_pars
                ),
            class = c("ForestTrainParam", "ForestTrain"))
}



#### New classifier ####
Param_classifier <- function(rastData, n_pts, Model, parallel=FALSE, progress = c("none","text")){

  thres <- as.numeric(Model$thres)
  matdata <- list(
                  as.matrix(rastData[[1]], wide=T),
                  as.matrix(rastData[[2]], wide=T),
                  as.matrix(rastData[[3]], wide=T)
                 )
  Gst <- Group_stats(matdata=matdata, n_pts=n_pts, fun=MultipleTesterParam2, progress = progress,
                     parallel = parallel, paropts=list(.packages = c("deforestable")), Model=Model)
  Gst_clfd <- classifier_green(data=Gst, thres=thres)

  cc <- max(Gst$cc)
  data_rest <- clsfd_matr_restor(c=cc, n_pts=n_pts, data=Gst_clfd)
  data_rest
}


