# library(doParallel)
# library(foreach)
# registerDoParallel()
library(deforestable)

test_that("ParamTrain works on data it was trained on", {
suppressWarnings({
  n_pts <- 20

  # Choosing folders with training data
  Forestdir <- system.file('extdata/Forest/', package = "deforestable")
  Nonforestdir <- system.file('extdata/Non-forest/', package = "deforestable")

  Model_P_tr <- train(model='fr_Param', Forestdir=Forestdir, Nonforestdir=Nonforestdir,
                      train_method='train', parallel=FALSE)

  # forest image 1
  test_image <- read_data_raster('_10_1_.jpeg', dir = Forestdir)
  res <- classify(data=test_image, Model=Model_P_tr,
                  n_pts=n_pts, parallel=FALSE, progress = 'none')
  expect_equal(sum(res), 400)

  # forest image 2
  test_image <- read_data_raster('_10_46_.jpeg', dir = Forestdir)
  res <- classify(data=test_image, Model=Model_P_tr,
                  n_pts=n_pts, parallel=FALSE, progress = 'none')
  expect_equal(sum(res), 400)

  # forest image 3
  test_image <- read_data_raster('_13_79_.jpeg', dir = Forestdir)
  res <- classify(data=test_image, Model=Model_P_tr,
                  n_pts=n_pts, parallel=FALSE, progress = 'none')
  expect_equal(sum(res), 400)

  # forest image 4
  test_image <- read_data_raster('_45_86_.jpeg', dir = Forestdir)
  res <- classify(data=test_image, Model=Model_P_tr,
                  n_pts=n_pts, parallel=FALSE, progress = 'none')
  expect_equal(sum(res), 400)

  # forest image 5
  test_image <- read_data_raster('_54_36_.jpeg', dir = Forestdir)
  res <- classify(data=test_image, Model=Model_P_tr,
                  n_pts=n_pts, parallel=FALSE, progress = 'none')
  expect_equal(sum(res), 400)

  # forest image 6
  test_image <- read_data_raster('_8_42_.jpeg', dir = Forestdir)
  res <- classify(data=test_image, Model=Model_P_tr,
                  n_pts=n_pts, parallel=FALSE, progress = 'none')
  expect_equal(sum(res), 400)

  ###### non-forest images

  # non-forest image 1
  test_image <- read_data_raster('_63_72_.jpeg', dir = Nonforestdir)
  res <- classify(data=test_image, Model=Model_P_tr,
                  n_pts=n_pts, parallel=FALSE, progress = 'none')
  expect_equal(sum(res), 0)

  # non-forest image 2
  test_image <- read_data_raster('_69_84_.jpeg', dir = Nonforestdir)
  res <- classify(data=test_image, Model=Model_P_tr,
                  n_pts=n_pts, parallel=FALSE, progress = 'none')
  expect_equal(sum(res), 0)

  # non-forest image 3
  test_image <- read_data_raster('_75_91_.jpeg', dir = Nonforestdir)
  res <- classify(data=test_image, Model=Model_P_tr,
                  n_pts=n_pts, parallel=FALSE, progress = 'none')
  expect_equal(sum(res), 0)

  # non-forest image 4
  test_image <- read_data_raster('_90_52_.jpeg', dir = Nonforestdir)
  res <- classify(data=test_image, Model=Model_P_tr,
                  n_pts=n_pts, parallel=FALSE, progress = 'none')
  expect_equal(sum(res), 0)

  # non-forest image 5
  test_image <- read_data_raster('_78_79_.jpeg', dir = Nonforestdir)
  res <- classify(data=test_image, Model=Model_P_tr,
                  n_pts=n_pts, parallel=FALSE, progress = 'none')
  expect_equal(sum(res), 0)

  # non-forest image 6
  test_image <- read_data_raster('_75_81_.jpeg', dir = Nonforestdir)
  res <- classify(data=test_image, Model=Model_P_tr,
                  n_pts=n_pts, parallel=FALSE, progress = 'none')
  expect_equal(sum(res), 0)
})
})













