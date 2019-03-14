#' Wrapper function for the autoencoder prediction + filtering step
#'
#' @inheritParams preprocessDat
#' @param use.pretrain Use a pretrained model or not
#' @param pretrained.weights.file If a pretrained model is used, provide the file storing the autoencoder model weights. It should have an extension of ".hdf5" and is the saved weights from the Python package \code{sctransfer}
#' @param save.ori Whether save the original.file.name to a new file
#' @param clearup.python Whether to clear up everything in the Python session after computation or not
#' @param ... more arguments passed to \code{autoFilterCV}
#' 
#' @return RDS file saved for the autoencoder prediction + filtering result
#' @export
computePrediction <- function(text.file.name, 
							  data.species = c("Human", "Mouse", "Others"), 
							  use.pretrain = F,
							  pretrained.weights.file = "",
							  model.species = c("Human", "Mouse", "Joint"),
							  model.nodes.ID = NULL, 
							  save.ori = T,
                clearup.python = T,
							  ...) {
	### inpute checking  ###
	format <- strsplit(text.file.name, '[.]')[[1]]
	format <- paste(".", format[length(format)], sep = "")
	if (format != ".txt" && format != ".csv" && format != ".rds")
		stop("Input file must be in .txt or .csv or .rds form", call.=FALSE)
	print(paste("Input file is:", text.file.name))
  if (use.pretrain)
    temp <- "Yes"
  else
    temp <- "No"
  print(paste("Use a pretrained model", temp))
	
	data.species <- match.arg(data.species, c("Human", "Mouse", "Others"))
  if (use.pretrain)
    print(paste("Data species is:", data.species)) 	

	if (use.pretrain) {
		if (data.species == "Others")
			stop("For pretrained model, the data.species can only be Human or Mouse")
		if (!file.exists(pretrained.weights.file))
			stop("Can not find the pretrained model. Please make sure that pretrained.weights.file exists")
		print(paste("Pretrained weights file is:", pretrained.weights.file))
		model.species <- match.arg(model.species, c("Human", "Mouse", "Joint"))
    if (model.species == "Joint")
      model.species <- data.species
		print(paste("Model species is:", model.species))

	}
	######

	### import Python module ###
	sctransfer <- reticulate::import("sctransfer", convert = F)
	main <- reticulate::import_main(convert = F)
	print("Python module sctransfer imported ...")
	######

	### preprocess data  ###
	if (use.pretrain && is.null(model.nodes.ID)) {
		if (model.species == "Human") {
	#		data(human_nodes_ID)
			model.nodes.ID <- human_nodes_ID
		} else {
	#		data(mouse_nodes_ID)
			model.nodes.ID <- mouse_nodes_ID
		}
	}
	preprocessDat(text.file.name, 
				  data.species = data.species, 
				  model.species = model.species, 
				  model.nodes.ID = model.nodes.ID, 
				  save.ori = save.ori)
	

	out_dir <- strsplit(text.file.name, split = "/")[[1]]
	out_dir <- paste(out_dir[-length(out_dir)], collapse = "/")
	if (out_dir == "")
		  out_dir <- "."
	print("Data preprocessed ...")
	######

	### run autoencoder ###
	if (use.pretrain) {
		x <- Matrix::readMM(gsub(format, ".mtx", text.file.name))
		x <- as.matrix(x)
		nonmissing_indicator <- read.table(gsub(format, "_nonmissing.txt", text.file.name))$V1
		batch.size <- as.integer(max(ncol(x) / 50, 32))
		used.time <- system.time(result <- autoFilterCV(x, 
												 sctransfer,
												 main,
												 pretrain_file = pretrained.weights.file,
												 nonmissing_indicator = nonmissing_indicator, 
												 model.species = model.species, 
												 out_dir = out_dir,
												 ...))
		print(paste("Autoencoder total computing time is:", used.time[3], "seconds"))

		idx <- nonmissing_indicator == 1
		print(paste("Number of predictive genes is", sum(result$err.const[idx] > result$err.autoencoder[idx])))

		data <- readRDS(gsub(format, ".rds", text.file.name))

		est.mu <- Matrix::rowMeans(Matrix::t(Matrix::t(data$mat) / Matrix::colSums(data$mat)) * 10000)
		est.mu <- est.mu %*% t(rep(1, ncol(data$mat)))

		rownames(result$x.autoencoder) <- model.nodes.ID

		if (data.species == model.species)
			ID.use <- data$rowdata$internal_ID
		else
			ID.use <- data$rowdata$other_species_internal_ID
		temp <- table(ID.use)
		id.not.unique <- names(temp[temp > 1])

		idx <- ID.use %in% model.nodes.ID[nonmissing_indicator == 1]
		est.mu[idx, ] <- result$x.autoencoder[ID.use[idx], ]
		names(result$err.autoencoder) <- names(result$err.const) <- rownames(result$x.autoencoder)
		err.autoencoder <- rep(NA, nrow(data$mat))
		err.const <- rep(NA, nrow(data$mat))
		err.autoencoder[idx] <- result$err.autoencoder[ID.use[idx]]
		err.const[idx] <- result$err.const[ID.use[idx]]
		file.remove(gsub(format, "_nonmissing.txt", text.file.name))
		file.remove(gsub(format, ".mtx", text.file.name))
		result$x.autoencoder <- est.mu
		result$err.autoencoder <- err.autoencoder
		result$err.const <- err.const
	} else {
		data <- readRDS(gsub(format, ".rds", text.file.name))

		x <- as.matrix(data$mat)
		batch.size <- as.integer(max(ncol(x) / 50, 32))
		used.time <- system.time(result <- autoFilterCV(x, 
												 sctransfer, 
												 main,
												 out_dir = out_dir, 
												 ...))

		print(paste("Autoencoder total computing time is:", used.time[3], "seconds"))
		est.mu <- result$x.autoencoder

		print(paste("Number of predictive genes is", sum(result$err.const > result$err.autoencoder)))

	}

	if (!use.pretrain || data.species == model.species) {
    temp.name <- gsub(format, "_prediction.rds", text.file.name)
		saveRDS(result, file = temp.name)
    print(paste("Predicted + filtered results saved as:", temp.name))
	} else {
    temp.name <- gsub(format, "_other_species_prediction.rds", text.file.name)
		saveRDS(result, file = temp.name)
    print(paste("Predicted + filtered results saved as:", temp.name))
	}

  reticulate::py_run_string("
import sys
sys.modules[__name__].__dict__.clear()")
  print("Python module cleared up.")
}








