#' Pre-process the data for SAVER-X 
#'
#'
#' @param out.dir a directory store all intemediate files. Can be generated automatically by this function is called by \code{computePrediction}
#' @param input.file.name Can be either .txt, .csv or .rds files that store the data matrix gene by cell, or can be NULL is \code{data.matrix} is provided. The rds file can either store the data as a regular matrix or a sparse matrix
#' @param data.matrix a matrix of UMI counts. Should be NULL if input.file.name is provided. The matrix is gene by cell with gene and cell names. Can be either a regular or sparse matrix.
#' @param data.species The species of the dataset
#' @param model.species The species of the pretrained model
#' @param model.nodes.ID The vector of node IDs of the pretrained model (only needed for the species of the data when the pre-trained model is joint). Set to NULL if running SAVER-X without pretraining.
#' 
#' @return Do not return anything but generate a series of intermediate files in the \code{out.dir}
#'
#' @export
preprocessDat <- function(out.dir,
                          input.file.name = NULL,
                          data.matrix = NULL,
                          data.species = c("Human", "Mouse", "Others"),
                          model.species = c("Human", "Mouse", "Joint"),
                          model.nodes.ID = NULL) {

  data.species <- match.arg(data.species, c("Human", "Mouse", "Others"))
  model.species <- match.arg(model.species, c("Human", "Mouse", "Joint"))

  if (is.null(input.file.name) && is.null(data.matrix))
    stop("Either an input data file or an input matrix should be provided!")

  if (!is.null(input.file.name) && !is.null(data.matrix))
    stop("Only either an input data file or an input matrix should be provided!")


#  task.id <- as.character(as.numeric(Sys.time()))  ## get a unique task ID to create new folder
#  dir.create(task.id)

  dir.create(out.dir, showWarnings = F)


  if (model.species == "Joint")
    model.species <- data.species

  if (!is.null(input.file.name)) {
    format <- strsplit(input.file.name, '[.]')[[1]]
    format <- paste(".", format[length(format)], sep = "")

    if (format == ".txt") {
      dat <- Matrix::Matrix(as.matrix(read.table(input.file.name)), sparse = T)
    } else if (format == ".csv") {
      dat <- Matrix::Matrix(as.matrix(read.csv(input.file.name, row.names = 1)), sparse = T)
    } else { 
      dat <- readRDS(input.file.name)
      if (is.list(dat))
        dat <- dat$mat
      dat <- Matrix::Matrix(dat, sparse = T)
    }
  } else {
    data.matrix <- Matrix::Matrix(data.matrix, sparse = T)
    dat <- data.matrix
    rm(data.matrix)
  }


  if (is.null(model.nodes.ID)) {
    dat <- list(mat = dat)
  #  temp.name <- gsub(format, "_temp.rds", text.file.name)
    temp.name <- paste0(out.dir, "/tmpdata.rds")
    saveRDS(dat, file = temp.name)
    print(paste("Processed file saved as:", temp.name))
  } else {
    if (is.null(rownames(dat)))
      stop("Need to provide gene names as row names of the data matrix")
    if (is.null(colnames(dat)))
      stop("Need to provide cell names as column names of the data matrix")

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
    if (sum(!is.na(result$internal.ID)) == 0)
      stop("No gene names can be mapped, please provide valid gene names!")

    print(paste(sum(!is.na(result$internal.ID)), "genes mapped out of", nrow(dat$mat)))

    temp.name <- paste0(out.dir, "/tmpdata.rds")
  #  temp.name <- gsub(format, "_temp.rds", text.file.name)
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

    if ("j" %in% methods::slotNames(mat1))
      mat <- Matrix::sparseMatrix(i = dat.node.idx[mat1@i+ 1], 
                          j = mat1@j, x = mat1@x, dims = c(length(model.nodes.ID), ncol(mat1)),
                          index1 = F)
    else
      mat <- Matrix::sparseMatrix(i = dat.node.idx[mat1@i+ 1], 
                          p = mat1@p, x = mat1@x, dims = c(length(model.nodes.ID), ncol(mat1)),
                          index1 = F)

    temp.name <- paste0(out.dir, "/tmpdata.mtx")
  #  temp.name <- gsub(format, ".mtx", text.file.name)
    Matrix::writeMM(mat, temp.name)
    print(paste("Reshaped file saved as:", temp.name))

    temp.name <- paste0(out.dir, "/tmpdata_nonmissing.txt")
#    temp.name <- gsub(format, "_nonmissing.txt", text.file.name)
    write.table(nonmissing.indicator, temp.name,
                row.names = F, col.names = F)
    print(paste("Nonmissing indicator saved as:", temp.name))

  }

}


