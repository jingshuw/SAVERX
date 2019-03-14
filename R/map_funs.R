#' map given mouse gene names to our internal unique 67618 mouse gene IDs
#'
#' @param gene.names a vector of given gene names that need to be mapped to
#'
#' @export
mapMouse <- function(gene.names) {

  gene.names <- as.character(gene.names)
  internal.ID <- rep(NA, length(gene.names))

  idx <- which(!(gene.names %in% rownames(all.dat.mouse)))
  
  if (length(idx) > 0) {
    internal.ID[-idx] <- gene.names[-idx]
    idx1 <- gene.names[idx] %in% rownames(alias.mouse)

    if (sum(idx1) > 0) {
      internal.ID[idx[idx1]] <- as.character(alias.mouse[gene.names[idx][idx1], "symbol"])
      idx <- idx[!idx1]
    }
  } else
    internal.ID <- gene.names

  if (length(idx) > 0) {
    idx1 <- gene.names[idx] %in% rownames(wiki.mouse)

    if (sum(idx1) > 0) {
      temp <- wiki.mouse[gene.names[idx][idx1], ]
      val <- temp[,"symbol"]
      val[val == ""] <- temp[val == "", "ensembl_gene_id"]

      internal.ID[idx[idx1]] <- val
      idx <- idx[!idx1]
    }
  }

  if (length(idx) > 0) {
    idx1 <- gene.names[idx] %in% rownames(moreID.mouse)

    if (sum(idx1) > 0) {
      internal.ID[idx[idx1]] <- as.character(moreID.mouse[gene.names[idx][idx1], "symbol"])
      idx <- idx[!idx1]
    }
  }


  symbol <- rep(NA, length(gene.names))
  ensembl_gene_id <- rep(NA, length(gene.names))
  if (length(idx) > 0) {
    symbol[-idx] <- as.character(all.dat.mouse[internal.ID[-idx], "symbol"])
    ensembl_gene_id[-idx] <- as.character(all.dat.mouse[internal.ID[-idx], "ensembl_gene_id"])
  } else {
    symbol <- as.character(all.dat.mouse[internal.ID, "symbol"])
    ensembl_gene_id <- as.character(all.dat.mouse[internal.ID, "ensembl_gene_id"])
  }

  return(list(internal.ID = internal.ID,
              symbol = symbol,
              ensembl_gene_id = ensembl_gene_id))

  
}

#' map given human gene names to our internal unique 64123 human gene IDs
#'
#' @param gene.names a vector of given gene names that need to be mapped to
#'
#' @export

mapHuman <- function(gene.names) {

  gene.names <- as.character(gene.names)
  internal.ID <- rep(NA, length(gene.names))

 
  idx <- which(!(toupper(gene.names) %in% names(human.upperID)))

  if (length(idx) > 0) {
    internal.ID[-idx] <- human.upperID[toupper(gene.names[-idx])]
    idx1 <- gene.names[idx] %in% rownames(alias.human)

    if (sum(idx1) > 0) {
      internal.ID[idx[idx1]] <- as.character(alias.human[gene.names[idx][idx1], "symbol"])
      idx <- idx[!idx1]
    }
  } else {
    internal.ID <- human.upperID[toupper(gene.names)]
  }

  if (length(idx) > 0) {
    idx1 <- toupper(gene.names[idx]) %in% rownames(prev.human)

    if (sum(idx1) > 0) {
      internal.ID[idx[idx1]] <- as.character(prev.human[toupper(gene.names[idx][idx1]), "symbol"])
      idx <- idx[!idx1]
    }
  }


  if (length(idx) > 0) {
    idx1 <- toupper(gene.names[idx]) %in% rownames(wiki.human)

    if (sum(idx1) > 0) {
      temp <- wiki.human[toupper(gene.names[idx][idx1]), ]
      val <- temp[,"symbol"]
      val[val == ""] <- temp[val == "", "ensembl_gene_id"]

      internal.ID[idx[idx1]] <- val
      idx <- idx[!idx1]
    }
  }

  if (length(idx) > 0) {
    idx1 <- gene.names[idx] %in% rownames(moreID.human)

    if (sum(idx1) > 0) {
      internal.ID[idx[idx1]] <- as.character(moreID.human[gene.names[idx][idx1], "symbol"])
      idx <- idx[!idx1]
    }
  }


  symbol <- rep(NA, length(gene.names))
  ensembl_gene_id <- rep(NA, length(gene.names))
  if (length(idx) > 0) {
    symbol[-idx] <- as.character(all.dat.human[internal.ID[-idx], "symbol"])
    ensembl_gene_id[-idx] <- as.character(all.dat.human[internal.ID[-idx], "ensembl_gene_id"])
  } else {
    symbol <- as.character(all.dat.human[internal.ID, "symbol"])
    ensembl_gene_id <- as.character(all.dat.human[internal.ID, "ensembl_gene_id"])
  }

  return(list(internal.ID = internal.ID,
              symbol = symbol,
              ensembl_gene_id = ensembl_gene_id))

  
}



