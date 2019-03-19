#' Cross-validation filtering of autoencoder
#'
#' Cross-validation is done to determine which genes can not be predicted well, by comparing the autoencoder predicted loss with the loss estimating the gene expression as a constant across cells
#'
#' @inheritParams autoencode
#' @param fold Number of total CV folds
#' @param samp Number of sampled folds taken to reduce CV time cost
#' @return a list of the filtered predicted data matrix and the CV error
#'
#' @export

autoFilterCV <- function(x,
                         python.module,
                         main,
                         pretrain_file = "",
                         nonmissing_indicator = 1,
                         n_human=21183L,
                         n_mouse=21122L,
                         shared_size=15494L,
                         model.species = NULL,
                         out_dir = ".",
                         batch_size = 32L,
                         write_output_to_tsv = F,
                         fold = 6, samp = 3, epsilon = 1e-10, seed = 1, ...) {

  set.seed(seed)

  n.cell <- ncol(x)
  n.gene <- nrow(x)

  idx.perm <- sample(1:n.cell, n.cell)

  n.test <- floor(n.cell/fold)
  
  err.autoencoder <- err.const <- rep(0, n.gene)
  ss <- sample(1:fold, samp)
  for (k in 1:samp) {
    print(paste("Cross-validation round:", k))
    i <- ss[k]
    train.idx <- idx.perm[-(((i-1)*n.test + 1):(i * n.test))]
    test.idx <- idx.perm[((i-1)*n.test + 1):(i * n.test)]

    x.test <- x[, test.idx]
    x.autoencoder <- autoencode(x[, train.idx],  
                                python.module,
                                main,
                                x.test,
                                pretrain_file,
                                nonmissing_indicator,
                                n_human,
                                n_mouse,
                                shared_size,
                                model.species,
                                out_dir,
                                batch_size, 
                                write_output_to_tsv,
                                verbose_sum = F, verbose_fit = 0L, ...)

    if (write_output_to_tsv) {
      x.autoencoder <- t(as.matrix(data.table::fread(paste0(out_dir, 
                                                          "/SAVERX_temp_pred_mean_norm.tsv"), header = F)))
    }
                 
    est.mu <- Matrix::rowMeans(Matrix::t(Matrix::t(x[, train.idx]) / Matrix::colSums(x[, train.idx])) * 10000)
    est.autoencoder <- Matrix::t(Matrix::t(x.autoencoder) * Matrix::colSums(x.test)) / 10000 
    rm(x.autoencoder)
    gc()
    est.const <- est.mu %*% Matrix::t(Matrix::colSums(x.test)) / 10000
    err1 <- -Matrix::rowMeans(x.test * log(est.autoencoder + epsilon) - est.autoencoder)
    err2 <- -Matrix::rowMeans(x.test * log(est.const + epsilon) - est.const)
    
    err.autoencoder <- err.autoencoder + err1
    err.const <- err.const + err2
    rm(x.test, est.mu, est.autoencoder, est.const, err1, err2)
    gc()
  }

  est.mu <- Matrix::rowMeans(Matrix::t(Matrix::t(x) / Matrix::colSums(x)) * 10000)
#  est.const <- est.mu %*% t(rep(1, n.cell))
  gnames <- rownames(x)
  cnames <- colnames(x)

  print("Final prediction round using all cells. See below the summary of the autoencoder model:")
  x.autoencoder <- autoencode(x,  
                              python.module,
                              main,
                              NULL,
                              pretrain_file,
                              nonmissing_indicator,
                              n_human,
                              n_mouse,
                              shared_size,
                              model.species,
                              out_dir,
                              batch_size, 
                              write_output_to_tsv, 
                              verbose_fit = 0L, ...)
  rm(x)
  gc()

  if (write_output_to_tsv) {
    x.autoencoder <- t(as.matrix(data.table::fread(paste0(out_dir, 
                                                        "/SAVERX_temp_mean_norm.tsv"), header = F)))
    rownames(x.autoencoder) <- gnames
    colnames(x.autoencoder) <- cnames
  }
  
  x.autoencoder[err.autoencoder - err.const > 0, ] <- est.mu[err.autoencoder - err.const > 0]


  return(list(x.autoencoder = x.autoencoder, 
              err.autoencoder = err.autoencoder, 
              err.const = err.const))
}
