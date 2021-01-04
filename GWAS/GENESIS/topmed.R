
.calcMAF <- function(gds, sample.id) {
  seqSetFilter(gds, sample.id=sample.id, verbose=FALSE)
  ref.freq <- alleleFrequency(gds)
  pmin(ref.freq, 1-ref.freq)
}


.calcMAC <- function(gds, sample.id) {
  seqSetFilter(gds, sample.id=sample.id, verbose=FALSE)
  ref.cnt <- alleleCount(gds)
  n.obs <- SeqVarTools:::.nSampObserved(gds)
  round(pmin(ref.cnt, 2*n.obs - ref.cnt))
}

filterByMAF <- function(gds, sample.id=NULL, maf.min=0, verbose=TRUE) {
  stopifnot(maf.min >= 0 & maf.min <= 0.5)
  if (maf.min == 0) return(invisible())
  if (sum(seqGetFilter(gds)$variant.sel) == 0) return(invisible())
  
  if (is.null(sample.id)) sample.id <- seqGetData(gds, "sample.id")
  maf <- .calcMAF(gds, sample.id)
  maf.filt <- maf >= maf.min
  if (verbose) message(paste("Running on", sum(maf.filt), "variants with MAF >=", maf.min))
  seqSetFilter(gds, variant.sel=maf.filt, action="intersect", verbose=verbose)
}

filterByMAC <- function(gds, sample.id=NULL, mac.min=1, verbose=TRUE) {
  stopifnot(mac.min >= 0)
  if (mac.min == 0) return(invisible())
  if (sum(seqGetFilter(gds)$variant.sel) == 0) return(invisible())
  
  if (is.null(sample.id)) sample.id <- seqGetData(gds, "sample.id")
  mac <- .calcMAC(gds, sample.id)
  maf.filt <- mac >= mac.min
  if (verbose) message(paste("Running on", sum(maf.filt), "variants with MAC >=", mac.min))
  seqSetFilter(gds, variant.sel=maf.filt, action="intersect", verbose=verbose)
}

checkSelectedVariants <- function(gds) {
  nvar <- sum(seqGetFilter(gds)$variant.sel)
  if (nvar == 0) {
    message("No variants selected. Exiting gracefully.")
    q(save="no", status=0)
  } else {
    message("Selected ", nvar, " variants.")
  }
}

filterByPass <- function(gds, verbose=TRUE) {
  filt <- seqGetData(gds, "annotation/filter")
  seqSetFilter(gds, variant.sel=(filt == "PASS"), action="intersect", verbose=verbose)
}

getobj <- function(Rdata) {
  objname <- load(Rdata)
  if (length(objname) > 1) {
    warning(paste("Multiple objects stored in file", Rdata,
                  "\nReturning only the first object"))
  }
  return(get(objname))
}

checkSelectedVariants <- function(gds) {
  nvar <- sum(seqGetFilter(gds)$variant.sel)
  if (nvar == 0) {
    message("No variants selected. Exiting gracefully.")
    q(save="no", status=0)
  } else {
    message("Selected ", nvar, " variants.")
  }
}

addMAC <- function(assoc, assoc_type) {
  mac <- function(x) {
    round(2 * x$n.obs * pmin(x$freq, 1-x$freq))
  }
  if (assoc_type == "single") {
    assoc$MAC <- mac(assoc)
  } else if (assoc_type %in% c("aggregate", "window")) {
    assoc$results$MAC <- sapply(assoc$variantInfo, function(x) sum(mac(x)))
  }
  assoc
}

countThreads <- function() {
  nSlots <- Sys.getenv("NSLOTS")
  nThreads <- ifelse(is.na(strtoi(nSlots) >= 1), 1, strtoi(nSlots))
  if (nThreads == 0) nThreads <- 1
  message(paste("Running with", nThreads,"thread(s)."))
  Sys.setenv(MKL_NUM_THREADS=nThreads)
  nThreads
}


##GRM functions

getGRM <- function(grm_file, sample.id=NULL) {
  files <- .splitFiles(grm_file)
  grm <- lapply(files, .readGRM, sample.id, matrix.name="grm")
  
  return(grm)
}

.splitFiles <- function(f) {
  strsplit(f, " ", fixed=TRUE)[[1]]
}

.readGRM <- function(f, sample.id, matrix.name="grm") {
  if (tools::file_ext(f) == "gds") {
    x <- openfn.gds(f)
    samp <- read.gdsn(index.gdsn(x, "sample.id"))
    if (is.null(sample.id)) sample.id <- samp
    sel <- samp %in% sample.id
    grm <- readex.gdsn(index.gdsn(x, matrix.name), sel=list(sel,sel))
    colnames(grm) <- rownames(grm) <- samp[sel]
    closefn.gds(x)
  } else {
    x <- getobj(f)
    if (matrix.name %in% names(x)) {
      colnames(x[[matrix.name]]) <- rownames(x[[matrix.name]]) <- x$sample.id
      if (!is.null(sample.id)) {
        keep <- x$sample.id %in% sample.id
        grm <- x[[matrix.name]][keep,keep]
      } else {
        grm <- x
      }
    } else {
      if (!is.null(sample.id)) {
        keep <- colnames(x) %in% sample.id
        grm <- x[keep, keep]
      } else {
        grm <- x
      }
    }
  }
  grm
}

