# Non-parametric cross-validation
#
# @examples
# library(deforeStable)
# n_pts <- 10
#
# forestdir <- system.file('extdata/Forest/', package = "deforeStable")
# Nonforestdir <- system.file('extdata/Non-forest/', package = "deforeStable")
#
# Best_model <- ParamCV(k_folds=3, n_pts=n_pts, forestdir=forestdir, Nonforestdir=Nonforestdir, parallel=FALSE)
# @export
ParamCV <- function(k_folds=NULL, n_pts=n_pts, forestdir=forestdir, Nonforestdir=Nonforestdir, parallel=FALSE)
{

    #### Obtaining k folds #######################
    trainPart <- createFolds(forestdir=forestdir, Nonforestdir=Nonforestdir, k = k_folds)


    kFolds_F <- trainPart$forest
    F_full <- unlist(kFolds_F)

    kFolds_NF <- trainPart$nonforest
    NF_full <- unlist(kFolds_NF)
    ##############################################

    Models <- list()
    test_ress <- list()
    acc <- vector(mode='numeric')

    for (ind in 1:length(kFolds_F)){

        # Defining the training sets
        train_F <- F_full[! F_full %in% kFolds_F[[ind]]]
        train_NF <- NF_full[! NF_full %in% kFolds_NF[[ind]]]

        # Training the model
        Models[[ind]] <- ParamTrain(forestdir = forestdir, Nonforestdir = Nonforestdir,
                                    Forest_list = train_F, Non_Forest_list = train_NF)

        # classifying test forest images
        f_npar_pos <- vector(mode='numeric')
        f_npar_sum <- vector(mode='numeric')

        for (ind_f in 1:length(kFolds_F[[ind]])) {
          data <- read_data_raster(filename = kFolds_F[[ind]][ind_f],
                                   dir=forestdir)
          dd <- Param_classifier(rastData=data, n_pts=n_pts, progress = "none",
                                 Model = Models[[ind]], parallel=parallel)
          f_npar_pos <- c(f_npar_pos, sum(dd)) # positive pixels tagged
          f_npar_sum <- c(f_npar_sum, dim(dd)[1] * dim(dd)[2]) # pixels total
        }

        # classifying test non-forest images
        nf_npar_pos <- vector(mode='numeric')
        nf_npar_sum <- vector(mode='numeric')

        for (ind_nf in 1:length(kFolds_NF[[ind]])) {
          data <- read_data_raster(filename = kFolds_NF[[ind]][ind_nf],
                                   dir=Nonforestdir)
          dd <- Param_classifier(rastData=data, n_pts=n_pts, progress = "none",
                                 Model = Models[[ind]], parallel=parallel)
          nf_npar_pos <- c(nf_npar_pos, sum(dd)) # positive pixels tagged
          nf_npar_sum <- c(nf_npar_sum, dim(dd)[1] * dim(dd)[2]) # pixels total
        }
        acc[ind] <- (sum(f_npar_pos) + (sum(nf_npar_sum) - sum(nf_npar_pos))) /
                    (sum(f_npar_sum) + sum(nf_npar_sum))

        test_ress[[ind]] <- list(acc=acc[ind],
                               tp = sum(f_npar_pos),
                               fp = sum(nf_npar_pos),
                               tn = sum(nf_npar_sum) - sum(nf_npar_pos),
                               fn = sum(f_npar_sum) - sum(f_npar_pos)
                              )
    }

    # Pick and return the best model
    model_N <- which(acc==max(acc)); model_N <- min(model_N)
    Out_model <- Models[[model_N]]
    ress <- test_ress[[model_N]]
    Out_object <- structure(list(call = Out_model$call,
                                 acc = Out_model$acc,
                                 tp = Out_model$tp,
                                 fp = Out_model$fp,
                                 tn = Out_model$tn,
                                 fn = Out_model$fn,
                                 thres = Out_model$thres,
                                 f_stab_params = Out_model$f_stab_params,
                                 forest_sigmas = Out_model$forest_sigmas,
                                 inv_f_sigmas = Out_model$inv_f_sigmas,
                                 f_Z = Out_model$f_Z,
                                 t_pars = Out_model$t_pars,
                                 test_ress = ress
                                ),
                            class = c("ForestCVParam", "ForestTrainParam", "ForestTrain"))

    Out_object
}




##################


