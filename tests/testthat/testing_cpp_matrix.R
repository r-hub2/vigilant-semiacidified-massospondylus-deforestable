library(deforestable)
library(StableEstim)


test_that("ComplexCF's are the same", {

####
theta <- c(1.3, 0.25, 0.3, 1)
t_par <- 1
r1 <- deforestable:::Compute_charact_funct(t=t_par,  alpha=theta[1], beta=theta[2],
                                           sigma=theta[3], mu=theta[4])
r2 <- StableEstim:::ComputeComplexCF(t_par, alpha=theta[1], beta=theta[2],
                                     gamma=theta[3], delta=theta[4], pm=1)
expect_equal(r1, r2)
expect_equal(deforestable:::ComplexCF_cpp(t_par, theta),
             ComplexCF(t_par, theta, pm=1))

####
theta <- c(1, 0.25, 0.3, 1)
t_par <- 1
r1 <- deforestable:::Compute_charact_funct(t=t_par,  alpha=theta[1], beta=theta[2],
                                           sigma=theta[3], mu=theta[4])
r2 <- StableEstim:::ComputeComplexCF(t_par, alpha=theta[1], beta=theta[2],
                                     gamma=theta[3], delta=theta[4], pm=1)
expect_equal(r1, r2)
expect_equal(deforestable:::ComplexCF_cpp(t_par, theta),
             ComplexCF(t_par, theta, pm=1))

####
theta <- c(1, 0, 0.3, 1)
t_par <- 1
r1 <- deforestable:::Compute_charact_funct(t=t_par,  alpha=theta[1], beta=theta[2],
                                           sigma=theta[3], mu=theta[4])
r2 <- StableEstim:::ComputeComplexCF(t_par, alpha=theta[1], beta=theta[2],
                                     gamma=theta[3], delta=theta[4], pm=1)
expect_equal(r1, r2)
expect_equal(deforestable:::ComplexCF_cpp(t_par, theta),
             ComplexCF(t_par, theta, pm=1))

####
theta <- c(2, 0.25, 0.3, 0)
t_par <- 1
r1 <- deforestable:::Compute_charact_funct(t=t_par,  alpha=theta[1], beta=theta[2],
                                           sigma=theta[3], mu=theta[4])
r2 <- StableEstim:::ComputeComplexCF(t_par, alpha=theta[1], beta=theta[2],
                                     gamma=theta[3], delta=theta[4], pm=1)
expect_equal(r1, r2)
expect_equal(deforestable:::ComplexCF_cpp(t_par, theta),
             ComplexCF(t_par, theta, pm=1))
})




test_that('Parametric matrices are the same',{

  ####
  theta <- c(1.3, 0.25, 0.3, 1)
  t_par <- 1
  expect_equal(deforestable:::stbl_param_covmtrx_cpp(t_par, theta),
               deforestable:::stbl_param_covmtrx(t_par, theta))

  ####
  theta <- c(1, 0.25, 0.3, 1)
  t_par <- 1
  expect_equal(deforestable:::stbl_param_covmtrx_cpp(t_par, theta),
               deforestable:::stbl_param_covmtrx(t_par, theta))

  ####
  theta <- c(1, 0, 0.3, 1)
  t_par <- 1
  expect_equal(deforestable:::stbl_param_covmtrx_cpp(t_par, theta),
               deforestable:::stbl_param_covmtrx(t_par, theta))

  ####
  theta <- c(2, 0.25, 0.3, 0)
  t_par <- 1
  expect_equal(deforestable:::stbl_param_covmtrx_cpp(t_par, theta),
               deforestable:::stbl_param_covmtrx(t_par, theta))
})



######################################

#############

