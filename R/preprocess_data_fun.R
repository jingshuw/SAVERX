#' Pre-process the data for SAVER-X 
#'
#' Output and intemediate files are stored at the same directory as text.file.name. To avoid potential file name conflicts, please make the folder only contain text.file.name. DO NOT run two SAVER-X tasks for the data file in the same folder.
#'
#' @param text.file.name Can be either .txt, .csv or .rds files that store the data matrix gene by cell. The rds file can either store the data as a regular matrix or a sparse matrix
#' @param data.species The species of the dataset
#' @param model.species The species of the pretrained model
#' @param model.nodes.ID The vector of node IDs of the pretrained model (only needed for the species of the data when the pre-trained model is joint). Set to NULL if running SAVER-X without pretraining.
#' @param save.ori Save a copy of the original \code{text.file.name} or not 
#' 
#' @return Do not return anything but generate a series of intermediate files in the same folder as \code{text.file.name}
#'
#' @export
preprocessDat <- function(text.file.name,
                          data.species = c("Human", "Mouse", "Others"),
                          model.species = c("Human", "Mouse", "Joint"),
                          model.nodes.ID = NULL,
						  save.ori = T) {

  data.species <- match.arg(data.species, c("Human", "Mouse", "Others"))
  model.species <- match.arg(model.species, c("Human", "Mouse", "Joint"))


  if (model.species == "Joint")
    model.species <- data.species

  format <- strsplit(text.file.name, '[.]')[[1]]
  format <- paste(".", format[length(format)], sep = "")

  

  if (format == ".txt") {
    dat <- Matrix::Matrix(as.matrix(read.table(text.file.name)), sparse = T)
  } else if (format == ".csv") {
    dat <- Matrix::Matrix(as.matrix(read.csv(text.file.name, row.names = 1)), sparse = T)
  } else { 
    dat <- readRDS(text.file.name)
  	if (save.ori) {
      temp.name <- gsub(".rds", "_original.rds", text.file.name)
		  system(paste("cp", text.file.name, temp.name))
      print(paste("Original file saved temporarily as:", temp.name))
    }
	  if (is.list(dat))
      dat <- dat$mat
    dat <- Matrix::Matrix(dat, sparse = T)
  }


  if (is.null(model.nodes.ID)) {
    dat <- list(mat = dat)
    temp.name <- gsub(format, ".rds", text.file.name)
    saveRDS(dat, file = temp.name)
    print(paste("Processed file saved as:", temp.name))
  } else {
    gene.names <- rownames(dat)
    rowdata <- data.frame(feature_symbol = gene.names)
    coldata <- data.frame(total_features = Matrix::colSums(dat != 0))
    rownames(coldata) <- colnames(dat)

    dat <- list(mat = dat, coldata= coldata, rowdata = rowdata)
    
    if (data.species == "Mouse") {
      result <- mapMouse(gene.names)
      dat$rowdata$mgi_symbol <- result$symbol
      dat$rowdata$other_species_internal_ID <- all.dat.mouse[result$internal.ID, "Human_ID"]
    } else { 
      result <- mapHuman(gene.names)
      dat$rowdata$hgnc_symbol <- result$symbol
      dat$rowdata$other_species_internal_ID <- all.dat.human[result$internal.ID, "Mouse_ID"]
    }

    dat$rowdata$ensembl_gene_id <- result$ensembl_gene_id
    dat$rowdata$internal_ID <- result$internal.ID

    print(paste(sum(!is.na(result$internal.ID)), "genes mapped out of", nrow(dat$mat)))

    temp.name <- gsub(format, ".rds", text.file.name)
    saveRDS(dat, temp.name)
    print(paste("Gene names mapped, resulting file saved as:", temp.name))

    node.idx <- 1:length(model.nodes.ID)
    names(node.idx) <- model.nodes.ID

    if (data.species == model.species)
      ID.use <- dat$rowdata$internal_ID
    else
      ID.use <- dat$rowdata$other_species_internal_ID


    temp <- table(ID.use)
    id.not.unique <- names(temp[temp > 1])

    idx <- (!ID.use %in% id.not.unique) & (ID.use %in% model.nodes.ID)
    mat1 <- dat$mat[idx, ]
    rownames(mat1) <- ID.use[idx]

    nonmissing.indicator <- as.numeric(model.nodes.ID %in% rownames(mat1))

    dat.node.idx <- node.idx[rownames(mat1)] - 1

    if ("j" %in% slotNames(mat1))
      mat <- Matrix::sparseMatrix(i = dat.node.idx[mat1@i+ 1], 
                          j = mat1@j, x = mat1@x, dims = c(length(model.nodes.ID), ncol(mat1)),
                          index1 = F)
    else
      mat <- Matrix::sparseMatrix(i = dat.node.idx[mat1@i+ 1], 
                          p = mat1@p, x = mat1@x, dims = c(length(model.nodes.ID), ncol(mat1)),
                          index1 = F)

    temp.name <- gsub(format, ".mtx", text.file.name)
    Matrix::writeMM(mat, temp.name)
    print(paste("Reshaped file saved as:", temp.name))

    temp.name <- gsub(format, "_nonmissing.txt", text.file.name)
    write.table(nonmissing.indicator, temp.name,
                row.names = F, col.names = F)
    print(paste("Nonmissing indicator saved as:", temp.name))

  }

}


