#' Funcion para leer los datos fenotipicos y genotipicos para GWASPoly
#' @Parametros_del_Modelo
#' 
#' @ploidy = Ploidy level of the organism
#' @pheno.file = Phenotypic file in local folder
#' @geno.file = Genotypic file in local folder
#' @format = Format for the marker data. Formats are the same as in GWASpoly read.GWASpoly() function
#' @n.traits = numer of traits to analyze in the phenofile
#' @delim = delimiter used for both the genofile and phenofile. Default is "," for the csv file
#' @snp_remover = list of names of the markers that we want to rtemove from the genofile 
#' @colIDy = column identifier for the name of genotype that links the genofile and the phenofile
#' @imp = impute function. Default is mode, but it can be mean as well


GWAS_CSD <- function (ploidy, pheno.file, geno.file, format, n.traits, delim = ",", snp_remover = NULL, colIDy = NULL, imp = "mode") 
{
  s = Sys.time()
  if (is.null(colIDy)) {
    stop("Proveer numero de columna donde esta el ID del genotipo que lo liga con la base genetica.")
  }
  pheno <- read.table(file = pheno.file, header = T, as.is = T, 
                      check.names = F, sep = delim, na.strings = c("",NA))
  gid.pheno <- unique(pheno[, 1])
  
  
  if (format == "ACTG") {
    format <- "ACGT"
  }
  if (!is.element(format, c("AB", "numeric", "ACGT"))) {
    stop("Invalid genotype format.")
  }
  bases <- c("A", "C", "G", "T")
  get.ref <- function(x, format) {
    if (format == "numeric") {
      ref.alt <- c(0, 1)
    }
    if (format == "AB") {
      ref.alt <- c("A", "B")
    }
    if (format == "ACGT") {
      y <- paste(na.omit(x), collapse = "")
      ans <- apply(array(bases), 1, function(z, y) {
        length(grep(z, y, fixed = T))
      }, y)
      if (sum(ans) > 2) {
        stop("Error in genotype matrix: More than 2 alleles")
      }
      if (sum(ans) == 2) {
        ref.alt <- bases[which(ans == 1)]
      }
      if (sum(ans) == 1) {
        ref.alt <- c(bases[which(ans == 1)], NA)
      }
    }
    return(ref.alt)
  }

  geno <- data.table::fread(file = geno.file, 
                            header = T, sep = delim, 
                            data.table = F,na.strings = c("",NA))
  
  if (is.null(snp_remover)) {
    geno <- geno
    cat("No se removieron marcadores\n")
  } else {
    geno <- subset(geno, !geno[,1] %in% snp_remover)
    cat(paste0("Se removieron ", length(snp_remover)," marcadores\n"))
  }
  
  map <- data.frame(Marker = geno[, 1], Chrom = factor(geno[,2], ordered = T), 
                    Position = geno[, 3], stringsAsFactors = F)
  
  alelos <- grep(paste(c("REF","Reference","Referencia"), 
                       collapse = "|"), colnames(geno), ignore.case = T)
  alelos1 <- grep(paste(c("ALT","Alternative","Alternativos"), 
                        collapse = "|"), colnames(geno), ignore.case = T)
  
  if (purrr::is_empty(alelos) | purrr::is_empty(alelos1)) {
    markers <- as.matrix(geno[, -c(1:3)])
    rownames(markers) <- geno[, 1]
    gid.geno <- colnames(geno)[-c(1:3)]
    
    tmp <- apply(markers, 1, get.ref, format)
    map$Ref <- tmp[1, ]
    map$Alt <- tmp[2, ]
    } else {
    markers <- as.matrix(geno[, -c(1:3,alelos,alelos1)])
    rownames(markers) <- geno[, 1]
    gid.geno <- colnames(geno)[-c(1:3,alelos,alelos1)]
    map$Ref = geno[,alelos]
    map$Alt = geno[,alelos1]
    
  }
  
  if (is.element(format, c("AB", "ACGT"))) {
    M <- apply(cbind(map$Ref, markers), 1, function(x) {
      y <- gregexpr(pattern = x[1], text = x[-1], fixed = T)
      ans <- as.integer(lapply(y, function(z) {
        ifelse(z[1] < 0, ploidy, ploidy - length(z))
      }))
      return(ans)
    }) }  else {
    M <- t(markers)
  }
  rownames(M) <- gid.geno
  stopifnot(na.omit(M <= ploidy & M >= 0))
  MAF <- apply(M, 2, function(x) {
    AF <- mean(x, na.rm = T)/ploidy
    MAF <- ifelse(AF > 0.5, 1 - AF, AF)
  })
  polymorphic <- which(MAF > 0)
  M <- M[, polymorphic]
  map <- map[polymorphic, ]
  map <- map[order(map$Chrom, map$Position), ]
  M <- M[, map$Marker]
  
  monomorphic <- apply(M, 2,function(x) length(unique(x)))
  mon <- which(monomorphic != 1)
  mon2 <- names(mon)
  M <- M[, mon2]
  
  map <- map[map$Marker %in% mon2, ]
  map <- map[order(map$Chrom, map$Position), ]
  M <- M[, map$Marker]
  
  m <- nrow(map)
  cat(paste("Number of polymorphic markers:", m, "\n"))
  
  impute.mode <- function(x) {
    ix <- which(is.na(x))
    if (length(ix) > 0) {
      x[ix] <- as.integer(names(which.max(table(x))))
    }
    return(x)
  }
  
  impute.mean <- function(x) {
    ix <- which(is.na(x))
    if (length(ix) > 0) {
      x[ix] <- round(mean(x, na.rm = T),0)
    }
    return(x)
  }
  
  missing <- which(is.na(M))
  if (length(missing)>0) {
    if (any(as.integer(M)!=as.numeric(M),na.rm=T)) {
      #fractional values present
      cat("Missing marker data imputed with population mean \n")
      M <- apply(M,2,impute.mean)
    } else {
      cat("Missing marker data imputed with population mode \n")
      M <- apply(M,2,impute.mode)
    }
  }
  
  gid <- intersect(gid.pheno, gid.geno)
  pheno <- pheno[is.element(pheno[, 1], gid), ]
  M <- M[gid, ]
  N <- length(gid)
  cat(paste("N =", N, "individuals with phenotypic and genotypic information \n"))
  n.fixed <- ncol(pheno) - n.traits - 1
  if (n.fixed > 0) {
    fixed <- data.frame(pheno[, (n.traits + 2):ncol(pheno)], 
                        stringsAsFactors = F)
    fixed.names <- colnames(pheno)[(n.traits + 2):ncol(pheno)]
    colnames(fixed) <- fixed.names
    pheno <- data.frame(pheno[, 1:(1 + n.traits)], stringsAsFactors = F)
    cat(paste("Detected following fixed effects:\n", paste(fixed.names, 
                                                           collapse = "\n"), "\n", sep = ""))
  }
  else {
    fixed <- data.frame(NULL)
  }
  traits <- colnames(pheno)[-1]
  cat(paste("Detected following traits:\n", paste(traits, 
                                                  collapse = "\n"), "\n", sep = ""))
  return(new("GWASpoly.data", map = map, pheno = pheno, fixed = fixed, 
             geno = M, ploidy = ploidy))
  
  e = Sys.time()
  cat("\n El procedimiento duro ", round(e - s, 2), units(e-s))
}
