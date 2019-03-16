#' Wrapper function for the empirical Bayes shrinkage step
#'
#' @inheritParams preprocessDat
#' @param ncores number of cores that can be used for the SAVER shrinkage
#' 
#' @return RDS saved for the final denoised results
#' @export
computeShrinkage <- function(text.file.name, ncores = 1) {

	### check input ###
	format <- strsplit(text.file.name, '[.]')[[1]]
	format <- paste(".", format[length(format)], sep = "")
	if (format != ".txt" && format != ".csv" && format != ".rds")
		stop("Input file must be in .txt or .csv or .rds form", call.=FALSE)

	out_dir <- strsplit(text.file.name, split = "/")[[1]]
	out_dir <- paste(out_dir[-length(out_dir)], collapse = "/")
	if (out_dir == "")
		out_dir <- "."
	######

	### preprocess data ###
	data <- readRDS(gsub(format, "_temp.rds", text.file.name))
	result <- readRDS(gsub(format, "_prediction.rds", text.file.name))
	est.mu <- result$x.autoencoder
	pred <- result$err.const > result$err.autoencoder
	names(pred) <- rownames(data$mat)
	rm(result)

	if (file.exists(gsub(format, "_other_species_prediction.rds", text.file.name))) {
		result.other <- readRDS(gsub(format, "_other_species_prediction.rds", text.file.name))

		est.mu[!is.na(data$rowdata$other_species_internal_ID), ] <- 
			result.other$x.autoencoder[!is.na(data$rowdata$other_species_internal_ID), ]

		pred[!is.na(data$rowdata$other_species_internal_ID)] <- 
			result.other$err.const > result.other$err.autoencoder

		rm(result.other)
	}
	######

	### run SAVER shrinkage ###
	used.time <- system.time(x.autoencoder.saver <- SAVER::saver(data$mat, mu = est.mu, 
																 ncores = ncores))
	print(paste("Empirical Bayes shrinkage total computing time is:", used.time[3], "seconds"))
	x.autoencoder.saver$est.before.shrinkage <- est.mu
	x.autoencoder.saver$predictable <- pred
#	mat <- x.autoencoder.saver$estimate
  temp.name <- gsub(format, "_denoised.rds", text.file.name)
	saveRDS(x.autoencoder.saver, temp.name)
  print(paste("Final denoised results saved as:", temp.name))

	try(file.remove(gsub(format, "_prediction.rds", text.file.name)), silent = T)
  if (file.exists(gsub(format, "_other_species_prediction.rds", text.file.name)))
	  file.remove(gsub(format, "_other_species_prediction.rds", text.file.name))
	try(file.remove(gsub(format, "_temp.rds", text.file.name)), silent = T)
	try(file.remove(paste0(out_dir, "/weights.hdf5")), silent = T)
  print("Intermediate files removed. Finished!!")
	######
}





