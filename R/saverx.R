#' the SAVERX function
#'
#' Output and intemediate files are stored at the same directory as the input data file \code{text.file.name}. To avoid potential file name conflicts, please make the folder only contain \code{text.file.name}. DO NOT run two SAVER-X tasks for the data file in the same folder.
#'
#' @inheritParams computeShrinkage
#' @inheritParams computePrediction
#' @param verbose Whether to show more autoencoder optimization progress or not
#' @return name of the final denoised RDS data file
#' @export
saverx <- function(input.file.name = NULL, 
                   data.matrix = NULL, 
                   data.species = c("Human", "Mouse", "Others"), 
                   use.pretrain = F,
                   pretrained.weights.file = "",
                   model.species = c("Human", "Mouse", "Joint"),
                   model.nodes.ID = NULL,
                   is.large.data = F,
                   ncores = 1, 
                   verbose = F, batch_size = NULL, 
                   clearup.python.session = T, ...) {


	if (use.pretrain) {
		data.species <- match.arg(data.species, c("Human", "Mouse", "Others"))
		model.species <- match.arg(model.species, c("Human", "Mouse", "Joint"))
    if (model.species == "Joint")
      model.species <- data.species
    if (use.pretrain && (data.species != model.species))
      clearup.python.session <- F
	}

  task.id <- as.character(as.numeric(Sys.time()))

  computePrediction(task.id, input.file.name, data.matrix,
                    data.species,
                    use.pretrain, pretrained.weights.file,
					  model.species, model.nodes.ID, verbose = verbose,
            is.large.data = is.large.data, batch_size = batch_size,
            clearup.python.session = clearup.python.session, ...)


  if (use.pretrain && (data.species != model.species)) {
    print("For cross-species training, compute another round with no pre-pretraining model ...")
    computePrediction(task.id, input.file.name, data.matrix, data.species, 
                      use.pretrain = F, save.ori = F, verbose = verbose, 
                      batch_size = batch_size, 
                      is.large.data = is.large.data, ...)
  }		


  print("Start the empirical Bayes shrinkage step using SAVER ...")
	dat.name <- computeShrinkage(task.id, ncores = ncores)

  return(dat.name)

}
