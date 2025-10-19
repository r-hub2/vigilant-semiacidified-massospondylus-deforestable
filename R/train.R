#' Train models for forest detection
#'
#' As input data, the function needs two folders- Nonforestdir with images of non-forest and Forestdir with ones of forest.
#' train() uses all images in both folders to train a model. Putting an image into an incorrect folder is equivalent to
#' tagging the image incorrectly.
#'
#' Currently, both fr_Non-Param and fr_Param use parameter n_pts only in the testing part of cross-validation, not during training.
#' Training is always done on whole original images in the training folders.
#'
#' @param model which model to train
#' @param Forestdir path to the directory with (only) forest images
#' @param Nonforestdir path to the directory with (only) non-forest images
#' @param train_method how to train the model: simple training, cross-validation.
#' @param n_pts matters only when train_method='cv'. Defines the size of the square sub-frames into which images would be split during cross-validation.
#' @param k_folds matters only when train_method='cv'. The number of folds in the k-fold cross-validation setup.
#' @param parallel matters only when train_method='cv'. Boolean. whether or not use a parallel setting during cross-validation
#'
#' @return object of class ForestTrain potentially with a sub-class. See \code{\link{Class_ForestTrain}}.
#'
#' @examples
#' library(deforestable)
#' n_pts <- 20
#'
#' # Choosing folders with training data
#' Forestdir <- system.file('extdata/Forest/', package = "deforestable")
#' Nonforestdir <- system.file('extdata/Non-forest/', package = "deforestable")
#'
#' k_folds=3;
#'
#' #### Read the target image ####
#' tg_dir <- system.file('extdata/', package = "deforestable")
#' test_image <- read_data_raster('smpl_1.jpeg', dir = tg_dir)
#'
#'
#' #### Models ####
#'
#'
#' # Simple training of the non-parametric model
#' Model_nonP_tr <- train(model='fr_Non-Param', Forestdir=Forestdir, Nonforestdir=Nonforestdir,
#'                        train_method='train', parallel=FALSE)
#'
#' res <- classify(data=test_image, Model=Model_nonP_tr,
#'                 n_pts=n_pts, parallel=FALSE, progress = 'text')
#'
#' tmp_d <- tempdir(); tmp_d
#' jpeg::writeJPEG(image=res, target=paste(tmp_d,'Model_nonP_tr.jpeg', sep='/'))
#'
#' \donttest{
#' # Cross-validation of the non-parametric model
#' Model_nonP_cv <- train(n_pts=n_pts, model='fr_Non-Param', Forestdir=Forestdir,
#'                        Nonforestdir=Nonforestdir, train_method='cv',
#'                        k_folds=k_folds, parallel=FALSE)
#'
#' res <- classify(data=test_image, Model=Model_nonP_cv,
#'                 n_pts=n_pts, parallel=FALSE, progress = 'text')
#'
#' tmp_d <- tempdir(); tmp_d
#' jpeg::writeJPEG(image=res, target=paste(tmp_d,'Model_nonP_cv.jpeg', sep='/'))
#'
#' }
#'
#' \donttest{
#' # Cross-validation of the parametric model
#' Model_P_cv <- train(n_pts=n_pts, model='fr_Param', Forestdir=Forestdir,
#'                     Nonforestdir=Nonforestdir, train_method='cv',
#'                     k_folds=k_folds, parallel=FALSE)
#'
#' res <- classify(data=test_image, Model=Model_P_cv,
#'                 n_pts=n_pts, parallel=FALSE, progress = 'text')
#'
#' tmp_d <- tempdir(); tmp_d
#' jpeg::writeJPEG(image=res, target=paste(tmp_d,'Model_P_cv.jpeg', sep='/'))
#'
#'
#' # Simple training of the parametric model
#' Model_P_tr <- train(model='fr_Param', Forestdir=Forestdir, Nonforestdir=Nonforestdir,
#'                     train_method='train', parallel=FALSE)
#' res <- classify(data=test_image, Model=Model_P_tr,
#'                 n_pts=n_pts, parallel=FALSE, progress = 'text')
#'
#' tmp_d <- tempdir(); tmp_d
#' jpeg::writeJPEG(image=res, target=paste(tmp_d,'Model_P_tr.jpeg', sep='/'))
#' }
#'
#' @export
train <- function(n_pts, model=c('fr_Non-Param', 'fr_Param'), Forestdir, Nonforestdir,
                  train_method=c('cv', 'train'), k_folds, parallel=FALSE)
{
    ####
    Forest_dir_ls <- list.files(path=Forestdir)
    Non_Forest_dir_ls <- list.files(path=Nonforestdir)

    if(length(Forest_dir_ls)==0) stop('There is no forest image to train on')
    ####

    if(model=='fr_Non-Param'){

        if(train_method=='cv') {

            if(length(Forest_dir_ls) < k_folds) stop('There are not enough forest images to split into folds')
            if(length(Non_Forest_dir_ls) < k_folds) stop('There are not enough non-forest images to split into folds')

            res <- NonParamCV(k_folds=k_folds, n_pts=n_pts, forestdir=Forestdir,
                              Nonforestdir=Nonforestdir, parallel=parallel)
        } else if(train_method=='train') {
            res <- NonParamTrain(forestdir=Forestdir, Nonforestdir=Nonforestdir,
                                 Forest_list=NULL, Non_Forest_list=NULL)
        }

    } else if(model=='fr_Param') {

        if(train_method=='train') {
            res <- ParamTrain(forestdir=Forestdir, Nonforestdir=Nonforestdir,
                              Forest_list=NULL, Non_Forest_list=NULL)

        } else if(train_method=='cv') {

            if(length(Forest_dir_ls) < k_folds) stop('There are not enough forest images to split into folds')
            if(length(Non_Forest_dir_ls) < k_folds) stop('There are not enough non-forest images to split into folds')

            res <- ParamCV(k_folds=k_folds, n_pts=n_pts,
                           forestdir=Forestdir, Nonforestdir=Nonforestdir,
                           parallel=parallel)
        }
    } else {stop('The specified model is not supported')}

    res
}




