#' Wrapper function for the empirical Bayes shrinkage step
#'
#' @inheritParams preprocessDat
#' @param ncores number of cores that can be used for the SAVER shrinkage
#' @param gene.block.size number of genes in one block if auto.split = T
#' @param cell.threshold the minimum threshold of the size of data for considering auto split
#' @param whether perform auto.split or not for the gateway RM-shared partition
#' 
#' @return RDS saved for the final denoised results
#' @export
computeShrinkage <- function(out.dir, ncores = 1, 
							 gene.block.size = 10000, cell.threshold = 20000, auto.split = T) {

	### check input ###
#	format <- strsplit(text.file.name, '[.]')[[1]]
#	format <- paste(".", format[length(format)], sep = "")
#	if (format != ".txt" && format != ".csv" && format != ".rds")
#		stop("Input file must be in .txt or .csv or .rds form", call.=FALSE)
#
#	out_dir <- strsplit(text.file.name, split = "/")[[1]]
#	out_dir <- paste(out_dir[-length(out_dir)], collapse = "/")
#	if (out_dir == "")
#		out_dir <- ".
	######

	### preprocess data ###
	data <- readRDS(paste0(out.dir, "/tmpdata.rds"))
	result <- readRDS(paste0(out.dir, "/prediction.rds"))
	est.mu <- result$x.autoencoder
	pred <- result$err.const > result$err.autoencoder
	names(pred) <- rownames(data$mat)
	rm(result)

	if (file.exists(paste0(out.dir, "/other_species_prediction.rds"))) {
		result.other <- readRDS(paste0(out.dir, "/other_species_prediction.rds"))

		est.mu[!is.na(data$rowdata$other_species_internal_ID), ] <- 
			result.other$x.autoencoder[!is.na(data$rowdata$other_species_internal_ID), ]

		pred[!is.na(data$rowdata$other_species_internal_ID)] <- 
			result.other$err.const > result.other$err.autoencoder

		rm(result.other)
	}
	######

	### run SAVER shrinkage ###
	ncells <- ncol(data$mat)
	ngenes <- nrow(data$mat)
    temp.name <- paste0(out.dir, "/denoised")
	if (ncores > 1 && ngenes > gene.block.size && ncells > cell.threshold && auto.split) {
		used.time <- system.time({	
			system(paste0('mv ', out.dir, '/prediction.rds ', out.dir, '/denoised_est_before_shrinkage.rds'))
			system(paste0('rm ', out.dir, '/prediction.rds'))
			System.sleep(60)
			#	saveRDS(est.mu, paste0(temp.name, "_est_before_shrinkage.rds"))
			sf <- Matrix::colSums(data$mat)
			sf <- sf / mean(sf)
			rd <- ceiling(ngenes / gene.block.size)
			x.autoencoder.saver <- SAVER::saver(data$mat[1:gene.block.size, ],
												 mu = est.mu[1:gene.block.size, ],
												 size.factor = sf, ncores = ncores)
			x.autoencoder.saver$predictable <- pred[1:gene.block.size]
			saveRDS(x.autoencoder.saver, paste0(temp.name, "_genes_1_", gene.block.size, ".rds"))
			x.autoencoder.saver$info <- NULL
			for (r in 2:rd) {
				i.start <- (r - 1) * gene.block.size + 1
			    i.end <- min(r * gene.block.size, ngenes)	
				x.autoencoder.saver <- SAVER::saver(data$mat[i.start:i.end, ],
									 mu = est.mu[i.start:i.end, ],
									 size.factor = sf, ncores = ncores)
				x.autoencoder.saver$predictable <- pred[i.start:i.end]
				saveRDS(x.autoencoder.saver, paste0(temp.name, "_genes_", i.start, "_", i.end, ".rds"))	
				print(paste("Final denoised results saved as:", paste0(temp.name, "_*.rds")))
			}
		}) 
	} else { 
		used.time <- system.time({
			x.autoencoder.saver <- SAVER::saver(data$mat, mu = est.mu, 
												ncores = ncores)

			x.autoencoder.saver$est.before.shrinkage <- est.mu
			x.autoencoder.saver$predictable <- pred
			saveRDS(x.autoencoder.saver, paste0(temp.name, ".rds"))
			print(paste("Final denoised results saved as:", paste0(temp.name, ".rds")))
		})
	}
	print(paste("Empirical Bayes shrinkage total computing time is:", used.time[3], "seconds"))

	tmp <- suppressWarnings(file.remove(paste0(out.dir, "/prediction.rds")))
	tmp <- suppressWarnings(file.remove(paste0(out.dir, "/other_species_prediction.rds")))
	tmp <- suppressWarnings(file.remove(paste0(out.dir, "/tmpdata.rds")))
	tmp <- suppressWarnings(file.remove(paste0(out.dir, "/weights.hdf5")))
  print("Intermediate files removed. Finished!!")
	######
  return(temp.name)
}





