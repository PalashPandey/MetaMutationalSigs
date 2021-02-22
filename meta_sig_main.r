# "=================================================================
# MetaMutationalSigs
# Author: Palash Pandey (pp535@drexel.edu)
# Desc:
#   extract - extract signatures by either automatic or semi-automatic way.
# 			Of note, when you use manual way, you need to run 2 times,
# 			firstly you should set --manual to get signature estimation results,
# 			and secondly you should set --manual --number N to get N signatures.
# Usage:
#   MetaMutationalSigs --inputdir=<input_directory> [--genome=<genome>]
# Options:
#   -h --help     Show help message.
#   --version     Show version.
#   -i <file>, --input <file>       input VCF directory path
#   -g <genome>, --genome <genome>  can be hg19, hg38 or mm10, [default: hg19].
# =================================================================" -> doc

# if (!suppressMessages(require("docopt"))) {
#   install.packages("docopt", repos = "https://cloud.r-project.org")
# }

# library("docopt")
# arguments <- docopt(doc, "MetaMutationalSigs v1.0\n")
# print(arguments)

# ## Stop error parsing
# if (!exists("arguments")) {
#   quit("no", status = -1)
# }

# message(
#   "
# =================================================================
#   __  __      _        __  __       _        _   _                   _  _____ _
#  |  \/  |    | |      |  \/  |     | |      | | (_)                 | |/ ____(_)
#  | \  / | ___| |_ __ _| \  / |_   _| |_ __ _| |_ _  ___  _ __   __ _| | (___  _  __ _ ___
#  | |\/| |/ _ \ __/ _` | |\/| | | | | __/ _` | __| |/ _ \| '_ \ / _` | |\___ \| |/ _` / __|
#  | |  | |  __/ || (_| | |  | | |_| | || (_| | |_| | (_) | | | | (_| | |____) | | (_| \__ \
#  |_|  |_|\___|\__\__,_|_|  |_|\__,_|\__\__,_|\__|_|\___/|_| |_|\__,_|_|_____/|_|\__, |___/
# 																				 __/ |
# 																				|___/
# Name   :       MetaMutationalSigs
# Link   :       https://github.com/PalashPandey/MetaMutationalSigs
# Doc    :       https://github.com/PalashPandey/MetaMutationalSigs
# ============================== START
# "
# )

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library(MutationalPatterns)
library(deconstructSigs)
library(sigfit)
library(sigminer)
data("cosmic_signatures_v2")
data("cosmic_signatures_v3")
load("data/signatures.exome.cosmic.v3.may2019.rda")
load("data/signatures.cosmic.rda")


# Parsing input
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("At least one argument must be supplied (input directory).",
       call. = FALSE)
} else if (length(args) == 1) {
  # default reference genome
  args[2] =  "hg19"
}

input <- args[1]
genome_build = args[2]

setwd(
  "C:\\Users\\pande\\OneDrive - Drexel University\\Documents\\Fall-2021\\Coop\\CGC\\SanjeeVCFFiles/PLOS_review_paper/metaSignatures/"
)

input = "test_data/maf/tcga_laml.maf.gz"
genome_build = "hg19"


print(input)
print(genome_build)

if (!is.null(input)) {
  if (dir.exists(input)) {
    fs <- list.files(input, pattern = "*.vcf", full.names = TRUE)
    if (length(fs) < 1) {
      message("When input is a directory, it should only contain VCF files!")
      quit("no", status = -1)
    }
    message("Try parsing VCF files...")
    obj <-
      sigminer::read_vcf(
        fs,
        genome_build = genome_build,
        keep_only_pass = FALSE,
        verbose = TRUE
      )
  } else {
    message("Try reading as a MAF file.")
    if (!is.data.frame(input)) {
      obj <-
        suppressMessages(data.table::fread(input, header = TRUE, data.table = FALSE))
    } else {
      obj <- input
    }
    if (!"Tumor_Sample_Barcode" %in% colnames(obj)) {
      message(
        "No 'Tumor_Sample_Barcode' column detected, try parsing as a component-by-sample matrix."
      )
      rownames(obj) <- obj[[1]]
      obj[[1]] <- NULL
      obj <- as.matrix(obj)
      message("Read as a matrix done.")
    } else {
      obj <- tryCatch(
        sigminer::read_maf(obj, verbose = TRUE),
        error = function(e) {
          message("Read input as a MAF file failed, try parsing as a component-by-sample matrix.")
          rownames(obj) <- obj[[1]]
          obj[[1]] <- NULL
          obj <- as.matrix(obj)
          message("Read as a matrix done.")
          return(obj)
        }
      )
      
    }
  }
}

data_matrix = sig_tally(obj)
write.csv(obj@data, "rectal_normal.maf")
data_matrix = data_matrix$nmf_matrix

data_matrix_id = sig_tally(obj, mode = c("ID"))
data_matrix_id = data_matrix_id$nmf_matrix

data_matrix_dbs = sig_tally(obj, mode = c("DBS"))
data_matrix_dbs = data_matrix_dbs$nmf_matrix

######### Installation and setup
FrobeniusNorm <- function(M, P, E) {
  sqrt(sum((M - P %*% E) ^ 2))
}
FrobeniusNorm_deconstructsigs <-
  function(i,
           input_matrices_list,
           exposures_list ,
           referenence_sigs) {
    M <- as.numeric(unlist(input_matrices_list[i]))
    P <- as.numeric(unlist(exposures_list[i]))
    E <- as.matrix(referenence_sigs)
    sqrt(sum((M - (P %*% E)) ^ 2))
  }
FrobeniusNorm_mutational_patterns <- function(M) {
  sqrt(sum((M) ^ 2))
}

whichSignaturesLocal = function(tumor.ref = NA,
                                sample.id,
                                signatures.ref = signatures.nature2013,
                                associated = c(),
                                signatures.limit = NA,
                                signature.cutoff = 0.06,
                                contexts.needed = FALSE,
                                tri.counts.method = "default") {
  if (class(tumor.ref) == 'matrix') {
    stop(paste(
      'Input tumor.ref needs to be a data frame or location of input text file',
      sep = ''
    ))
  }
  
  if (ncol(signatures.ref) == 78 & tri.counts.method != 'default') {
    warning('Using default normalization for DBS signatures')
    tri.counts.method <- 'default'
  }
  tumor     <- tumor.ref
  if (missing(sample.id) && nrow(tumor) == 1) {
    sample.id = rownames(tumor)[1]
  }
  # Take patient id given
  tumor <- as.matrix(tumor)
  if (!sample.id %in% rownames(tumor)) {
    stop(paste(sample.id, " not found in rownames of tumor.ref", sep = ''))
  }
  tumor <- subset(tumor, rownames(tumor) == sample.id)
  # if(round(rowSums(tumor), digits = 1) != 1){
  #   stop(paste('Sample: ', sample.id, ' is not normalized\n', 'Consider using "contexts.needed = TRUE"', sep = ' '))
  # }
  
  #Read in Stratton signatures file
  
  if (exists("signatures.ref", mode = "list")) {
    signatures   <- signatures.ref
  } else {
    if (file.exists(signatures.ref)) {
      signatures <-
        utils::read.csv(
          signatures.ref,
          sep = "\t",
          header = TRUE,
          as.is = TRUE,
          check.names = FALSE
        )
    } else {
      print("signatures.ref is neither a file nor a loaded data frame")
    }
  }
  
  signatures    <- as.matrix(signatures)
  original.sigs <- signatures
  
  # Check column names are formatted correctly
  
  if (length(colnames(tumor)[colnames(tumor) %in% colnames(signatures)]) < length(colnames(signatures))) {
    # colnames(tumor) <- changeColumnNames(colnames(tumor))
    if (length(colnames(tumor)[colnames(tumor) %in% colnames(signatures)]) < length(colnames(signatures))) {
      stop("Check column names on input file")
    }
  }
  
  # Ensure that columns in tumor match the order of those in signatures
  tumor <- tumor[, colnames(signatures), drop = FALSE]
  
  #Take a subset of the signatures
  if (!is.null(associated)) {
    signatures <-
      signatures[rownames(signatures) %in% associated, , drop = FALSE]
  }
  
  if (is.na(signatures.limit)) {
    signatures.limit <- nrow(signatures)
  }
  
  # Remove signatures from possibilities if they have a "strong" peak not seen in the tumor sample
  zero.contexts   <- colnames(tumor)[tumor < 0.01]
  corr.sigs       <-
    which(signatures[, zero.contexts, drop = FALSE] >= 0.2, arr.ind = T)
  signatures      <-
    signatures[which(!rownames(signatures) %in% rownames(corr.sigs)), , drop = FALSE]
  #print(paste(rownames(corr.sigs), " not considered in the analysis.", sep = ""))
  
  #Set the weights matrix to 0
  weights         <-
    matrix(
      0,
      nrow = nrow(tumor),
      ncol = nrow(signatures),
      dimnames = list(rownames(tumor), rownames(signatures))
    )
  
  seed            <- findSeed(tumor, signatures)
  weights[seed]   <- 1
  w               <- weights * 10
  
  error_diff      <- Inf
  error_threshold <- 1e-3
  
  num <- 0
  while (error_diff > error_threshold) {
    num        <- num + 1
    #print(num)
    error_pre  <- getError(tumor, signatures, w)
    if (error_pre == 0) {
      break
    }
    #print(w)
    w          <-
      updateW_GR(tumor, signatures, w, signatures.limit = signatures.limit)
    error_post <- getError(tumor, signatures, w)
    #print(paste("old_error = ", error_pre, sep = ""))
    #print(paste("new_error = ", error_post, sep = ""))
    #print(w)
    error_diff <- (error_pre - error_post) / error_pre
  }
  
  weights <- w / sum(w)
  unknown <- 0
  
  ## filtering on a given threshold value (0.06 default)
  weights[weights < signature.cutoff] <- 0
  unknown <- 1 - sum(weights)
  
  product <- weights %*% signatures
  diff    <- tumor - product
  
  x       <-
    matrix(
      data = 0,
      nrow = 1,
      ncol = nrow(original.sigs),
      dimnames = list(rownames(weights), rownames(original.sigs))
    )
  x       <- data.frame(x)
  x[colnames(weights)] <- weights
  weights <- x
  
  out        <- list(weights, tumor, product, diff, unknown)
  names(out) <- c("weights", "tumor", "product", "diff", "unknown")
  return(out)
  
}

#########################################
# Mutational Patterns
#########################################
#' Run mutationalPatterns refitting pipeline for COSMIC 30 and SBS signatures
#' @param ref_genome Refernce genome version. Needs to be in the BSgenome library.
#' @param out_dir The output directory
#' @return Mutational catalog i.e matrix of counts of 96 mutational contexts
#' @example sample_mut_mat = runMuationalPatterns(input_dir = "input_vcfs", ref_genome= "BSgenome.Hsapiens.UCSC.hg19", out_dir = "mutationalSigsPlots")
runMuationalPatterns = function(ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                                out_dir,
                                mut_mat,
                                mut_mat_id,
                                mut_mat_dbs) {
  mut_mat <- data_matrix
  mut_mat <- t(mut_mat)
  ref_genome <- ref_genome
  # COSMIC 2019 SBS signatures
  signatures = get_known_signatures(source = "COSMIC")
  
  # Strict signature fitting iteratively reduces the number of signatures used for refitting by removing the signature with least contribution
  strict_refit <-
    fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.004)
  fit_res_strict <- strict_refit$fit_res
  mut_patterns_strict_errors <-
    as.matrix.data.frame(t(mut_mat - fit_res_strict$reconstructed))
  mut_patterns_strict_errors <-
    apply(mut_patterns_strict_errors  ,
          1,
          FrobeniusNorm_mutational_patterns)
  setwd(out_dir)
  dir.create("./mutational_patterns_results")
  write.csv(
    fit_res_strict$contribution,
    "./mutational_patterns_results/strict_sample_exposures.csv"
  )
  write.csv(
    mut_patterns_strict_errors,
    row.names = colnames(mut_mat),
    "./mutational_patterns_results/strict_sample_errors.csv"
  )
  # This is traditional decomposition
  fit_res <- fit_to_signatures(mut_mat, signatures)
  mut_patterns_errors <-
    as.matrix.data.frame(t(mut_mat - fit_res$reconstructed))
  mut_patterns_errors <-
    apply(mut_patterns_errors  , 1, FrobeniusNorm_mutational_patterns)
  write.csv(fit_res$contribution,
            "./mutational_patterns_results/sample_exposures.csv")
  write.csv(mut_patterns_errors,
            row.names = colnames(mut_mat),
            "./mutational_patterns_results/sample_errors.csv")
  
  # COSMIC 30 signatures
  signatures_legacy = t(cosmic_signatures_v2)
  # Strict signature fitting iteratively reduces the number of signatures used for refitting by removing the signature with least contribution
  strict_refit_legacy <-
    fit_to_signatures_strict(mut_mat, signatures_legacy, max_delta = 0.004)
  fit_res_strict_legacy <- strict_refit_legacy$fit_res
  mut_patterns_strict_legacy_errors <-
    as.matrix.data.frame(t(mut_mat - fit_res_strict_legacy$reconstructed))
  mut_patterns_strict_legacy_errors <-
    apply(mut_patterns_strict_legacy_errors  ,
          1,
          FrobeniusNorm_mutational_patterns)
  write.csv(
    fit_res_strict_legacy$contribution,
    "./mutational_patterns_results/legacy_strict_sample_exposures.csv"
  )
  write.csv(
    mut_patterns_strict_legacy_errors,
    row.names = colnames(mut_mat),
    "./mutational_patterns_results/legacy_strict_sample_errors.csv"
  )
  # This is traditional decomposition
  fit_res_legacy <- fit_to_signatures(mut_mat, signatures_legacy)
  mut_patterns_legacy_errors <-
    as.matrix.data.frame(t(mut_mat - fit_res_strict$reconstructed))
  mut_patterns_legacy_errors <-
    apply(mut_patterns_legacy_errors  ,
          1,
          FrobeniusNorm_mutational_patterns)
  write.csv(
    fit_res_legacy$contribution,
    "./mutational_patterns_results/legacy_sample_exposures.csv"
  )
  write.csv(
    mut_patterns_legacy_errors,
    row.names = colnames(mut_mat),
    "./mutational_patterns_results/legacy_sample_errors.csv"
  )
  
  #  ID signatures
  mut_mat_id = t(data_matrix_id)
  signatures_id = get_known_signatures(source = "COSMIC", muttype = "indel")
  strict_refit <- fit_to_signatures(mut_mat_id, signatures_id)
  mut_patterns_id_errors <-
    as.matrix.data.frame(t(mut_mat_id - strict_refit$reconstructed))
  mut_patterns_id_errors <-
    apply(mut_patterns_id_errors  ,
          1,
          FrobeniusNorm_mutational_patterns)
  write.csv(strict_refit$contribution,
            "./mutational_patterns_results/id_sample_exposures.csv")
  write.csv(
    mut_patterns_id_errors,
    row.names = colnames(mut_mat_id),
    "./mutational_patterns_results/id_sample_errors.csv"
  )
  
  
  # DBS signatures
  mut_mat_dbs = t(data_matrix_dbs)
  signatures_dbs = get_known_signatures(source = "COSMIC", muttype = "dbs")
  strict_refit <- fit_to_signatures(mut_mat_dbs, signatures_dbs)
  mut_patterns_dbs_errors <-
    as.matrix.data.frame(t(mut_mat_dbs - strict_refit$reconstructed))
  mut_patterns_dbs_errors <-
    apply(mut_patterns_dbs_errors  ,
          1,
          FrobeniusNorm_mutational_patterns)
  write.csv(
    strict_refit$contribution,
    "./mutational_patterns_results/dbs_sample_exposures.csv"
  )
  write.csv(
    mut_patterns_dbs_errors,
    row.names = colnames(mut_mat_dbs),
    "./mutational_patterns_results/dbs_sample_errors.csv"
  )
  
  return(mut_mat)
}

#########################################
# DeconstructSigs
#########################################
#' Run DeconstructSigs refitting pipeline for COSMIC 30 and SBS signatures
#' @param ref_genome Refernce genome version. Needs to be in the BSgenome library.
#' @param out_dir The output directory
#' @example runDeconstructSigs(input_dir = "input_vcfs", ref_genome= "BSgenome.Hsapiens.UCSC.hg19", out_dir = "DeconstructSigsPlots")
runDeconstructSigs = function(ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                              out_dir,
                              mut_mat,
                              mut_mat_id,
                              mut_mat_dbs) {
  mut_mat <- data_matrix
  v2_input_matrices_list <- mut_mat
  v2_input_matrices_list <-
    apply(v2_input_matrices_list, 1, as.data.frame)
  v2_input_matrices_list <- lapply(v2_input_matrices_list, t)
  v2_input_matrices_list <-
    lapply(v2_input_matrices_list, as.data.frame)
  v2_signatures_list <-
    lapply(
      v2_input_matrices_list,
      whichSignatures,
      contexts.needed = TRUE,
      signatures.ref = signatures.cosmic
    )
  v2_exposures_list <-
    lapply(v2_signatures_list, function(x)
      x$weights)
  v2_diff_list <- lapply(v2_signatures_list, function(x)
    x$diff)
  v2_error_list <-
    lapply(seq(length(v2_diff_list))  , function(i)
      FrobeniusNorm_deconstructsigs(
        i,
        v2_input_matrices_list,
        v2_exposures_list,
        signatures.cosmic
      ))
  #### SBS
  #   v3_input_matrices_list <- lapply(snp_files_list, vcf.to.sigs.input)
  v3_input_matrices_list <- v2_input_matrices_list
  
  v3_signatures_list <-
    lapply(
      v3_input_matrices_list,
      whichSignatures,
      contexts.needed = TRUE,
      signatures.ref = signatures.exome.cosmic.v3.may2019
    )
  v3_exposures_list <-
    lapply(v3_signatures_list, function(x)
      x$weights)
  v3_diff_list <- lapply(v3_signatures_list, function(x)
    x$diff)
  v3_error_list <-
    lapply(seq(length(v3_diff_list))  , function(i)
      FrobeniusNorm_deconstructsigs(
        i,
        v3_input_matrices_list,
        v3_exposures_list,
        signatures.exome.cosmic.v3.may2019
      ))
  
  #  ID
  indel_signatures = get_known_signatures(source = "COSMIC", muttype = "indel")
  mut_mat_id = data_matrix_id
  indel_mutational_profile <- data.frame(mut_mat_id)
  indel_mutational_profile <-
    indel_mutational_profile %>% rename_at(colnames(indel_mutational_profile) ,
                                           ~ colnames(as.data.frame(t(indel_signatures))))
  indel_input_matrices_list <-
    lapply(1:dim(indel_mutational_profile)[1], function(x)
      as.data.frame(indel_mutational_profile[x[1], ]))
  indel_signatures_list <-
    lapply(
      indel_input_matrices_list,
      whichSignaturesLocal,
      contexts.needed = TRUE,
      signatures.ref = as.data.frame(t(indel_signatures))
    )
  indel_exposures_list <-
    lapply(indel_signatures_list, function(x)
      x$weights)
  indel_diff_list <-
    lapply(indel_signatures_list, function(x)
      x$diff)
  indel_error_list <-
    lapply(seq(length(indel_diff_list))  , function(i)
      FrobeniusNorm_deconstructsigs(
        i,
        indel_input_indel_matrices_list,
        indel_exposures_list,
        as.data.frame(t(indel_signatures))
      ))
  
  #  DBS
  dbs_signatures = get_known_signatures(source = "COSMIC", muttype = "dbs")
  mut_mat_dbs = data_matrix_dbs
  dbs_mutational_profile <- data.frame(mut_mat_dbs)
  dbs_mutational_profile <-
    dbs_mutational_profile %>% rename_at(colnames(dbs_mutational_profile) ,
                                         ~ colnames(as.data.frame(t(dbs_signatures))))
  dbs_input_matrices_list <-
    lapply(1:dim(dbs_mutational_profile)[1], function(x)
      as.data.frame(dbs_mutational_profile[x[1], ]))
  dbs_signatures_list <-
    lapply(
      dbs_input_matrices_list,
      whichSignaturesLocal,
      contexts.needed = TRUE,
      signatures.ref = as.data.frame(t(dbs_signatures))
    )
  dbs_exposures_list <-
    lapply(dbs_signatures_list, function(x)
      x$weights)
  dbs_diff_list <- lapply(dbs_signatures_list, function(x)
    x$diff)
  dbs_error_list <-
    lapply(seq(length(dbs_diff_list))  , function(i)
      FrobeniusNorm_deconstructsigs(
        i,
        dbs_input_matrices_list,
        dbs_exposures_list,
        as.data.frame(t(dbs_signatures))
      ))
  
  
  setwd(out_dir)
  dir.create("./deconstructsigs_results")
  write.csv(
    unlist(v2_error_list),
    row.names = rownames(data_matrix),
    "./deconstructsigs_results/legacy_sample_errors.csv"
  )
  write.csv(
    do.call(rbind, lapply(v2_exposures_list, data.frame)),
    row.names = rownames(data_matrix),
    "./deconstructsigs_results/legacy_sample_exposures.csv"
  )
  write.csv(
    unlist(v3_error_list),
    row.names = rownames(data_matrix),
    "./deconstructsigs_results/sbs_sample_errors.csv"
  )
  write.csv(
    do.call(rbind, lapply(v3_exposures_list, data.frame)),
    row.names = rownames(data_matrix),
    "./deconstructsigs_results/sbs_sample_exposures.csv"
  )
  
  write.csv(
    unlist(indel_error_list),
    row.names = rownames(mut_mat_id),
    "./deconstructsigs_results/indel_sample_errors.csv"
  )
  write.csv(
    do.call(rbind, lapply(indel_exposures_list, data.frame)),
    "./deconstructsigs_results/indel_sample_exposures.csv"
  )
  
  write.csv(
    unlist(dbs_error_list),
    row.names = rownames(mut_mat_dbs),
    "./deconstructsigs_results/dbs_sample_errors.csv"
  )
  write.csv(
    do.call(rbind, lapply(dbs_exposures_list, data.frame)),
    "./deconstructsigs_results/dbs_sample_exposures.csv"
  )
}

#########################################
# Sigfit
#########################################
#' Run Sigfit refitting pipeline for COSMIC 30 and SBS signatures
#' @param ref_genome Refernce genome version. Needs to be in the BSgenome library.
#' @param out_dir The output directory
#' @param  mut_mat Mutational catalog i.e matrix of counts of 96 mutational contexts
#' @example runSigfit(input_dir = "input_vcfs", ref_genome= "BSgenome.Hsapiens.UCSC.hg19", out_dir = "SigfitPlots")
runSigfit = function(ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                     out_dir,
                     mut_mat,
                     mut_mat_id,
                     mut_mat_dbs) {
  mut_mat = data_matrix
  samples_fit <-
    fit_signatures(
      counts = mut_mat,
      signatures = cosmic_signatures_v2,
      iter = 3000,
      warmup = 1000,
      chains = 1,
      seed = 1756
    )
  sample_exposures <- retrieve_pars(samples_fit, "exposures")
  catalogue_matrix <- as.data.frame(mut_mat)
  sig_matrix <- t(as.matrix(samples_fit[["data"]][["signatures"]]))
  expo_mat <- t(as.matrix(sample_exposures$mean))
  errors_per_sample <-
    sapply(seq(ncol(expo_mat)), function(i)
      FrobeniusNorm(catalogue_matrix[i,], sig_matrix, expo_mat[, i]))
  rmse <- sqrt(mean(errors_per_sample ^ 2))
  errors_per_sample <- as.data.frame(errors_per_sample)
  row.names(errors_per_sample) <- row.names(sample_exposures$mean)
  #### SBS
  samples_fit_sbs <-
    fit_signatures(
      counts = mut_mat,
      signatures = cosmic_signatures_v3,
      iter = 3000,
      warmup = 1000,
      chains = 1,
      seed = 1756
    )
  sample_exposures_sbs <-
    retrieve_pars(samples_fit_sbs, "exposures")
  catalogue_matrix_sbs <- t(as.data.frame(mut_mat))
  sig_matrix_sbs <-
    t(as.matrix(samples_fit_sbs[["data"]][["signatures"]]))
  expo_mat_sbs <- t(as.matrix(sample_exposures_sbs$mean))
  errors_per_sample_sbs <-
    sapply(seq(ncol(expo_mat_sbs)), function(i)
      FrobeniusNorm(catalogue_matrix_sbs[i,], sig_matrix_sbs, expo_mat_sbs[, i]))
  rmse_sbs <- sqrt(mean(errors_per_sample_sbs ^ 2))
  errors_per_sample_sbs <- as.data.frame(errors_per_sample_sbs)
  row.names(errors_per_sample_sbs) <-
    row.names(sample_exposures_sbs$mean)
  # Indels
  mut_mat_id = data_matrix_id
  indel_mutational_profile = mut_mat_id
  indel_signatures = get_known_signatures(source = "COSMIC", muttype = "indel")
  indel_samples_fit <-
    fit_signatures(
      counts = indel_mutational_profile ,
      signatures =  as.data.frame(t(indel_signatures)),
      iter = 3000,
      warmup = 1000,
      chains = 1,
      seed = 1756
    )
  indel_sample_exposures <-
    retrieve_pars(indel_samples_fit, "exposures")
  sig_matrix <-
    t(as.matrix(indel_samples_fit[["data"]][["signatures"]]))
  expo_mat <- t(as.matrix(indel_sample_exposures$mean))
  errors_per_sample <-
    sapply(seq(ncol(expo_mat)), function(i)
      FrobeniusNorm(indel_mutational_profile[i,], sig_matrix, expo_mat[, i]))
  errors_per_sample <- as.data.frame(errors_per_sample)
  row.names(errors_per_sample) <-
    row.names(indel_sample_exposures$mean)
  
  # DBS
  mut_mat_dbs = data_matrix_dbs
  dbs_mutational_profile = mut_mat_dbs
  dbs_signatures = get_known_signatures(source = "COSMIC", muttype = "dbs")
  dbs_samples_fit <-
    fit_signatures(
      counts = dbs_mutational_profile ,
      signatures =  as.data.frame(t(dbs_signatures)),
      iter = 3000,
      warmup = 1000,
      chains = 1,
      seed = 1756
    )
  dbs_sample_exposures <-
    retrieve_pars(dbs_samples_fit, "exposures")
  sig_matrix <-
    t(as.matrix(dbs_samples_fit[["data"]][["signatures"]]))
  expo_mat <- t(as.matrix(dbs_sample_exposures$mean))
  errors_per_sample_dbs <-
    sapply(seq(ncol(expo_mat)), function(i)
      FrobeniusNorm(dbs_mutational_profile[i,], sig_matrix, expo_mat[, i]))
  errors_per_sample_dbs <- as.data.frame(errors_per_sample_dbs)
  row.names(errors_per_sample_dbs) <-
    row.names(dbs_sample_exposures$mean)
  
  setwd(out_dir)
  dir.create("./sigfit_results")
  write.csv(errors_per_sample,
            "./sigfit_results/sample_errors_indel.csv")
  write.csv(indel_sample_exposures$mean,
            "./sigfit_results/sample_exposures_indel.csv")
  write.csv(errors_per_sample_dbs,
            "./sigfit_results/sample_errors_dbs.csv")
  write.csv(dbs_sample_exposures$mean,
            "./sigfit_results/sample_exposures_dbs.csv")
  write.csv(errors_per_sample_sbs,
            "./sigfit_results/sample_errors_sbs.csv")
  write.csv(sample_exposures_sbs$mean,
            "./sigfit_results/sample_exposures_sbs.csv")
  write.csv(errors_per_sample,
            "./sigfit_results/legacy_sample_errors.csv")
  write.csv(sample_exposures$mean, "./sigfit_results/legacy_sample_exposures.csv")
  
  
  
}
#########################################
# Sigflow
#########################################
#' Run Sigminer/Sigflow refitting pipeline for COSMIC 30 and SBS signatures
#' @param input_dir The directory where sample VCF files are
#' @param ref_genome Refernce genome version. Hg19/Hg38
#' @param out_dir The output directory
#' @example runSigflow(input_dir = "input_vcfs", ref_genome= "BSgenome.Hsapiens.UCSC.hg19", out_dir = "SigflowPlots")
runSigflow = function(ref_genome = "hg19",
                      out_dir,
                      mut_mat,
                      mut_mat_id,
                      mut_mat_dbs) {
  data_matrix = mut_mat
  indel_data_matrix = mut_mat_id
  dbs_data_matrix = mut_mat_dbs
  fit_legacy <-  sig_fit(
    catalogue_matrix = t(data_matrix),
    sig_index = "ALL",
    sig_db = "legacy",
    mode = "SBS",
    method = "SA",
    return_class = "data.table",
    return_error = TRUE
  )
  setwd(out_dir)
  output_fit(fit_legacy, result_dir = "sigflow" , mut_type = "legacy")
  
  
  fit_sbs <-  sig_fit(
    catalogue_matrix = t(data_matrix),
    sig_index = "ALL",
    sig_db = "SBS",
    mode = "SBS",
    method = "SA",
    return_class = "data.table",
    return_error = TRUE
  )
  output_fit(fit_sbs, result_dir =  "sigflow" , mut_type = "SBS")
  
  fit_sbs <-  sig_fit(
    catalogue_matrix = t(indel_data_matrix),
    sig_index = "ALL",
    sig_db = "ID",
    mode = "ID",
    method = "SA",
    return_class = "data.table",
    return_error = TRUE
  )
  output_fit(fit_sbs, result_dir =  "sigflow" , mut_type = "ID")
  
  fit_sbs <-  sig_fit(
    catalogue_matrix = t(dbs_data_matrix),
    sig_index = "ALL",
    sig_db = "DBS",
    mode = "DBS",
    method = "SA",
    return_class = "data.table",
    return_error = TRUE
  )
  output_fit(fit_sbs, result_dir =  "sigflow" , mut_type = "DBS")
}
#########################################
# Run commands
#########################################

output_dir = "C:\\Users\\pande\\OneDrive - Drexel University\\Documents\\Fall-2021\\Coop\\CGC\\SanjeeVCFFiles\\PLOS_review_paper\\meta_package_data"
output_dir = args[3]

print(output_dir)
sample_mut_mat = runMuationalPatterns(mut_mat = data_matrix, mut_mat_id = data_matrix_id, mut_mat_dbs = data_matrix_dbs,  out_dir = output_dir)
runDeconstructSigs(mut_mat = data_matrix, mut_mat_id = data_matrix_id, mut_mat_dbs = data_matrix_dbs, out_dir = output_dir)
runSigfit(mut_mat = data_matrix, mut_mat_id = data_matrix_id, mut_mat_dbs = data_matrix_dbs, out_dir = output_dir)
runSigflow(mut_mat = data_matrix, mut_mat_id = data_matrix_id, mut_mat_dbs = data_matrix_dbs, out_dir = output_dir)

# system(paste("python error_pie_heatmap.py " , output_dir))
