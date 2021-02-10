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

setwd("C:\\Users\\pande\\OneDrive - Drexel University\\Documents\\Fall-2021\\Coop\\CGC\\SanjeeVCFFiles\\PLOS_review_paper/metaSignatures/")

input = "test_data"
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
    obj <- sigminer::read_vcf(fs, genome_build = genome_build, keep_only_pass = FALSE, verbose = TRUE)
  } else {
    message("Try reading as a MAF file.")
    if (!is.data.frame(input)) {
      obj <- suppressMessages(data.table::fread(input, header = TRUE, data.table = FALSE))
    } else {
      obj <- input
    }
    if (!"Tumor_Sample_Barcode" %in% colnames(obj)) {
      message("No 'Tumor_Sample_Barcode' column detected, try parsing as a component-by-sample matrix.")
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
data_matrix = data_matrix$nmf_matrix
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
								mut_mat) {
  mut_mat <- t(mut_mat)
  ref_genome <- ref_genome
  # COSMIC 2019 SBS signatures
  signatures = get_known_signatures(source = "COSMIC")
  
  # Strict signature fitting iteratively reduces the number of signatures used for refitting by removing the signature with least contribution
  strict_refit <- fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.004)
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
    "./mutational_patterns_results/legacy_sample_errors.csv"
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
							  mut_mat) {
	v2_input_matrices_list <- mut_mat
	v2_input_matrices_list <- apply(v2_input_matrices_list, 1, as.data.frame)
	v2_input_matrices_list <- lapply(v2_input_matrices_list, t)
	v2_input_matrices_list <- lapply(v2_input_matrices_list, as.data.frame)
  v2_signatures_list <-  lapply(v2_input_matrices_list,whichSignatures,contexts.needed = TRUE,signatures.ref = signatures.cosmic)
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
                     mut_mat) {
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
      FrobeniusNorm(catalogue_matrix[i, ], sig_matrix, expo_mat[, i]))
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
      FrobeniusNorm(catalogue_matrix_sbs[i, ], sig_matrix_sbs, expo_mat_sbs[, i]))
  rmse_sbs <- sqrt(mean(errors_per_sample_sbs ^ 2))
  errors_per_sample_sbs <- as.data.frame(errors_per_sample_sbs)
  row.names(errors_per_sample_sbs) <-
    row.names(sample_exposures_sbs$mean)
  
  setwd(out_dir)
  dir.create("./sigfit_results")
  write.csv(errors_per_sample, "./sigfit_results/sample_errors_legacy.csv")
  write.csv(sample_exposures, "./sigfit_results/sample_exposures_legacy.csv")
  write.csv(errors_per_sample_sbs,
            "./sigfit_results/sample_errors_sbs.csv")
  write.csv(sample_exposures_sbs,
            "./sigfit_results/sample_exposures_sbs.csv")
  
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
					  mut_mat) {
  data_matrix = mut_mat 
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
}
#########################################
# Run commands
#########################################

output_dir = "C:\\Users\\pande\\OneDrive - Drexel University\\Documents\\Fall-2021\\Coop\\CGC\\SanjeeVCFFiles\\PLOS_review_paper\\meta_package_data"
output_dir = args[3]

print(output_dir)
sample_mut_mat = runMuationalPatterns(mut_mat = data_matrix, out_dir = output_dir)
runSigfit(mut_mat = data_matrix, out_dir = output_dir)
runDeconstructSigs(mut_mat = data_matrix, out_dir = output_dir)
runSigflow(mut_mat = data_matrix, out_dir = output_dir)

# system(paste("python error_pie_heatmap.py " , output_dir))
