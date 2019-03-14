#' R interface for the autoencode python function
#'
#' @param x Target sparse data matrix of gene by cell. When pretraining is used, the genes should be the same as the nodes used in the pretrained model. If a node gene is missing is the target dataset, set all values of that gene as 0 in \code{x} and indicate that using \code{nonmissing_indicator}
#' @param python.module The python module for the Python package \code{sctransfer}
#' @param main A Python main module
#' @param nonmissing_indicator A single value 1 or a vector of 0 and 1s to indicate which nodes are missing in the target dataset. Set to 1 for no pretraining.
#' @param test.x Data matrix to evaluate the test error
#' @param pretrain_file The pretrained weights file ended with '.hdf5'
#' @param  model.species Should be either 'Human' or 'Mouse' when pretraining is used
#' @param ... Extra parameters passed to Python module \code{sctransfer} function \code{api} (if no pretraining) or function \code{api_pretrain} (with pretraining).
#'
#' @return A data matrix for the Python autoencoder result 
#'
#' @export

#require(reticulate)

autoencode <- function(x, 
                       python.module,
                       main,
                       test.x = NULL,  
                       pretrain_file = "", 
                       nonmissing_indicator = 1, 
                       n_human=21183L,
                       n_mouse=21122L,
                       shared_size=15494L,
                       model.species = NULL,
                       out_dir = ".",
                       batch_size = 32L,
                       ...) {

  if (pretrain_file == "")
    pretrain <- F
  else
    pretrain <- T

  if (pretrain)
    api <- python.module$api_pretrain
  else
    api <- python.module$api
  

  gnames <- rownames(x)
  cnames <- colnames(x)
  x <- Matrix::Matrix(x, sparse = T)
  mtx_file <- paste0(out_dir, "/SAVERX_temp.mtx")
  Matrix::writeMM(x, file = mtx_file)
 # x <- api$anndata$AnnData(t(x))
 # main$x <- x
  nonmissing_indicator <- api$np$asarray(nonmissing_indicator) 

  if (!is.null(test.x)) {
    gnames <- rownames(test.x)
    cnames <- colnames(test.x)
    test.x <- api$anndata$AnnData(t(as.matrix(test.x)))
    main$test_x <- test.x
  }

  if (!pretrain)
    main$result <- api$autoencode(mtx_file = mtx_file,
                                  pred_adata = test.x,
                                  nonmissing_indicator = nonmissing_indicator,                      
                                  out_dir = out_dir,
                                  batch_size = batch_size, 
                                  ...)
  else
    main$result <- api$autoencode(n_inoutnodes_human=n_human,
                                  n_inoutnodes_mouse=n_mouse,
                                  shared_size=shared_size,
                                 # adata = x,
                                  mtx_file = mtx_file,
                                  pred_adata = test.x, 
                                  species=model.species,
                                  nonmissing_indicator = nonmissing_indicator, 
                                  initial_file = pretrain_file, 
                                  out_dir = out_dir,
                                  batch_size = batch_size,
                                  ...)

  x.autoencoder <- t(reticulate::py_to_r(main$result$obsm[['X_dca']]))
  colnames(x.autoencoder) <- cnames
  rownames(x.autoencoder) <- gnames
  reticulate::py_run_string("
del result
import gc
gc.collect()")
  return(x.autoencoder)
}
