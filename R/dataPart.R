#' Data Partitioning
#'
#' As input data, the functions need two folders- Nonforestdir with images of non-forest and forestdir with ones of forest.
#' createDataPartition() splits data into training and testing partitions while keeping the relative sample size of the classes
#' the same as in the original data. createFolds() splits the data into k folds for cross-validation.
#'
#' @param forestdir path to the directory with (only) forest images
#' @param Nonforestdir path to the directory with (only) non-forest images
#' @param times the number of data partitions to make
#' @param p the percentage of data to set aside for training
#' @param k the number of folds to split the data into
#'
#' @return createDataPartition returns a list of data partitions. Each partition consists of 4 sets- forest training,
#' non-forest training, forest test and non-forest test set. createFolds returns lists $forest and $nonforest with k folds in 
#' each of them.
#'
#' @examples
#'
#' library(deforestable)
#' forestdir <- system.file('extdata/Forest/', package = "deforestable")
#' Nonforestdir <- system.file('extdata/Non-forest/', package = "deforestable")
#'
#' trainPart <- createDataPartition(forestdir=forestdir, Nonforestdir=Nonforestdir, p = .7, times = 1)
#'
#' folds <- createFolds(forestdir, Nonforestdir, k = 10)
#'
#' @export
createDataPartition <- function (forestdir, Nonforestdir, times = 1, p = 0.5){

  Forest_list <- list.files(path=forestdir)
  Non_Forest_list <- list.files(path=Nonforestdir)

  # Combine image names into matrix with Forest and Not Forest labels
  dataset <- rbind(cbind(Im = Forest_list, label = rep("Forest", length(Forest_list))),
                   cbind(Im = Non_Forest_list, label = rep("Not Forest", length(Non_Forest_list))))

  output <- vector(mode = "list", times)

  if(nrow(dataset) < 2){
    stop("Input has too few data points")
  }

  for(i in 1:times) {
    y <- dataset[sample(1:nrow(dataset)), ]
    yF <- y[y[,2]=="Forest",]
    yNF <- y[y[,2]=="Not Forest",]

    curr <- vector(mode = "list", 0)

    curr$ForTrain <- yF[1:ceiling(nrow(yF)*p),1]
    curr$NForTrain <- yNF[1:ceiling(nrow(yNF)*p),1]

    curr$ForTest <- yF[-(1:ceiling(nrow(yF)*p)),1]
    curr$NForTest <- yNF[-(1:ceiling(nrow(yNF)*p)),1]

    output[[i]] <- curr
  }
  output
}


#' @describeIn createDataPartition Split data into folds
#' @export
createFolds <- function(forestdir, Nonforestdir, k = 5){


  Forest_list <- list.files(path=forestdir)
  Non_Forest_list <- list.files(path=Nonforestdir)

  # Combine image names into matrix with Forest and Not Forest labels
  dataset <- rbind(cbind(Im = Forest_list, label = rep("Forest", length(Forest_list))),
                   cbind(Im = Non_Forest_list, label = rep("Not Forest", length(Non_Forest_list))))


  y <- dataset[sample(1:nrow(dataset)), ]
  yF <- y[y[,2]=="Forest",]
  yNF <- y[y[,2]=="Not Forest",]

  if(length(yF[,1]) < k | length(yNF[,1]) < k){
    stop("Less observations than folds in atleast one of the classes, reduce folds")
  }


  forestfolds <- split(yF[,1], cut(seq(1, nrow(yF)), breaks=k, labels=FALSE))
  nonforestfolds <- split(yNF[,1], cut(seq(1, nrow(yNF)), breaks=k, labels=FALSE))

  out <- list(forest = forestfolds, nonforest = nonforestfolds)

}

