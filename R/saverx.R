#' the SAVERX function
#'
#' Output and intemediate files are stored at the same directory as the input data file \code{text.file.name}. To avoid potential file name conflicts, please make the folder only contain \code{text.file.name}. DO NOT run two SAVER-X tasks for the data file in the same folder.
#'
#' @inheritParams computeShrinkage
#' @inheritParams computePrediction
#' @param verbose Whether to show more autoencoder optimization progress or not
#' @return the final denoised RDS data file saved in the same directory as the input data
#' @export
saverx <- function(text.file.name, 
				   data.species = c("Human", "Mouse", "Others"), 
				   use.pretrain = F,
				   pretrained.weights.file = "",
				   model.species = c("Human", "Mouse", "Joint"),
				   model.nodes.ID = NULL,
				   ncores = 1, 
           verbose = F, ...) {

	computePrediction(text.file.name, data.species,
					  use.pretrain, pretrained.weights.file,
					  model.species, model.nodes.ID, verbose = verbose, ...)

	if (use.pretrain) {
		data.species <- match.arg(data.species, c("Human", "Mouse"))
		model.species <- match.arg(model.species, c("Human", "Mouse", "Joint"))
    if (model.species == "Joint")
      model.species <- data.species
		if (data.species != model.species) {
			print("For cross-species training, compute another round with no pre-pretraining model ...")
			computePrediction(text.file.name, data.species, 
							  use.pretrain = F, save.ori = F, verbose = verbose, ...)
		}		
	}

  print("Start the empirical Bayes shrinkage step using SAVER ...")
	computeShrinkage(text.file.name, ncores = ncores)

}
