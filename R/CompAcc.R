CompAccNP <- function(ds_f, ds_nf){

  mi <- min(ds_nf)
  ma <- max(ds_f)

  if (mi>ma) {
    thres <- (mi+ma)/2
    err <- 0
    acc <- 1
    tp <- length(ds_f)
    fp <- 0
    tn <- length(ds_nf)
    fn <- 0

  } else {
    # how many points there are in the intersection
    fs <- ds_f[ds_f > mi]
    nfs<- ds_nf[ds_nf < ma]

    # sort these points and find thresholds
    vect <- sort(c(fs,nfs))
    out <- (vect[2:length(vect)] + vect[1:(length(vect)-1)])/2

    n_err <- vector(mode='numeric')
    n_err_f_vec <- vector(mode='numeric')
    n_err_nf_vec <- vector(mode='numeric')
    for (ind in seq(from=1, by=1, along.with=(out))) {
      n_err_f <- length(fs[fs > out[ind]])
      n_err_nf <- length(nfs[nfs < out[ind]])
      n_err <- c(n_err, n_err_f + n_err_nf)
      n_err_f_vec <- c(n_err_f_vec, n_err_f)
      n_err_nf_vec <- c(n_err_nf_vec, n_err_nf)
    }

    best_thres_ind <- which(n_err==min(n_err))
    best_thres_ind <- floor(mean(best_thres_ind))

    thres <- out[best_thres_ind]
    err <- (n_err[best_thres_ind] / (length(ds_f) + length(ds_nf)))[1]
    acc <- 1-err
    tp <- length(ds_f) - n_err_f_vec[best_thres_ind]
    fp <- n_err_nf_vec[best_thres_ind]
    tn <- length(ds_nf) - n_err_nf_vec[best_thres_ind]
    fn <- n_err_f_vec[best_thres_ind]
  }

  list(thres=thres, acc=acc, tp = tp, fp = fp, tn = tn, fn = fn)
}
