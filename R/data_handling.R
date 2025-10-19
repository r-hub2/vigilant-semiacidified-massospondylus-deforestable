
#' Import a jpeg image
#'
#' All these functions are made to read jpeg images, the difference is in the class of objects they return
#'
#' @return read_data returns a 3-column data.frame with pixels in rows and
#' red, green, blue intensities in columns. read_data_matrix reads jpeg images and returns 3 matrices for
#' each of red, green and blue colors. read_data_raster imports jpeg as a raster object
#' \code{\link[terra]{rast}}.
#'
#' @param filename name of the jpeg file to import
#' @param dir the directory where the image is located
#' @examples
#'
#' dir <- system.file('extdata/Forest/', package = "deforestable")
#'
#' dd <- read_data(filename='_6_33_.jpeg', dir=dir)
#' hist(dd[,1])
#'
#' @export
read_data <- function(filename, dir){

  fil <- paste(dir, filename, sep = "/")
  obj <- jpeg::readJPEG(source=fil)

  red <- as.vector(obj[,,1])
  green <- as.vector(obj[,,2])
  blue <- as.vector(obj[,,3])

  res <- data.frame(red, green, blue)
  return(res)
}


#' Import a jpeg image
#'
#' @examples
#'
#' dir <- system.file('extdata/Forest/', package = "deforestable")
#'
#' dd <- read_data_matrix(filename='_6_33_.jpeg', dir=dir)
#'
#' @describeIn read_data returns three matrices
#' @export
read_data_matrix <- function(filename, dir){

  fil <- paste(dir, filename, sep = "/")
  obj <- jpeg::readJPEG(source=fil, native = FALSE)

  red <- obj[,,1]
  green <- obj[,,2]
  blue <- obj[,,3]


  return(list(red, green, blue))
}


#' Import jpeg as a raster object
#'
#'
#' @examples
#'
#' dir <- system.file('extdata/Forest/', package = "deforestable")
#'
#' dd<-read_data_raster(filename='_8_46_.jpeg', dir=dir)
#'
#' @describeIn read_data returns a SpatRaster object
#' @export
read_data_raster <- function(filename, dir){

  fil <- paste(dir, filename, sep = "/")
  obj <- jpeg::readJPEG(source=fil)
  res <- terra::rast(obj)
  return(res)
}








