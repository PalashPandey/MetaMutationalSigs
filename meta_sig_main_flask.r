#!/usr/bin/env Rscript
# setwd("../")
library(MutationalPatterns)
library(deconstructSigs)
library(sigfit)
library(sigminer)
library(dplyr)
library(ggpubr)
library(tidyverse)
library("devtools")

######### Installation and setup
FrobeniusNorm <- function(M, P, E) {
  sqrt(sum((M - P %*% E) ^ 2))
}
FrobeniusNorm_deconstructsigs <- function(i, input_matrices_list, exposures_list , referenence_sigs) {
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
    error_pre  <- getError(tumor, signatures, w)
    if (error_pre == 0) {
      break
    }
    w          <-
      updateW_GR(tumor, signatures, w, signatures.limit = signatures.limit)
    error_post <- getError(tumor, signatures, w)
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
# Legacy
#########################################
runLegacy = function(ref_genome = "hg19",
                      mut_mat_input, mutationalPatterns = TRUE, sigflow = TRUE , sigfit = TRUE , deconstructSigs = TRUE) {
  mut_mat <- mut_mat_input
  mut_mat <- t(mut_mat)

  ####### MutationalPatteTRUE
  if (mutationalPatterns){
    signatures_legacy = as.matrix(legacy_signatures)
    # Strict signature fitting iteratively reduces the number of signatures used for refitting by removing the signature with least contribution
    strict_refit_legacy <-    fit_to_signatures_strict(mut_mat, signatures_legacy, max_delta = 0.004)
    fit_res_strict_legacy <- strict_refit_legacy$fit_res
    mut_patterns_strict_legacy_errors <- as.matrix.data.frame(t(mut_mat - fit_res_strict_legacy$reconstructed))
    mut_patterns_strict_legacy_errors <-apply(mut_patterns_strict_legacy_errors  , 1, FrobeniusNorm_mutational_patterns)
    dir.create("./mutational_patterns_results")
    write.csv(t(fit_res_strict_legacy$contribution),"./mutational_patterns_results/legacy_strict_sample_exposures.csv" )
    write.csv(mut_patterns_strict_legacy_errors,row.names = colnames(mut_mat),"./mutational_patterns_results/legacy_strict_sample_errors.csv")
    # This is traditional decomposition
    fit_res_legacy <- fit_to_signatures(mut_mat, signatures_legacy)
    mut_patterns_legacy_errors <-   as.matrix.data.frame(t(mut_mat - fit_res_legacy$reconstructed))
    mut_patterns_legacy_errors <- apply(mut_patterns_legacy_errors  ,1,FrobeniusNorm_mutational_patterns)
    write.csv(t(fit_res_legacy$contribution),"./mutational_patterns_results/legacy_sample_exposures.csv")
    write.csv(mut_patterns_legacy_errors,row.names = colnames(mut_mat),"./mutational_patterns_results/legacy_sample_errors.csv"  )
    mut_legacy_df = as_tibble(read.csv("mutational_patterns_results/legacy_strict_sample_errors.csv") )%>% 
      rename( sample = X,  error = x)%>%
      mutate(toolname = "mutationalPatterns")
    final_legacy_df = mut_legacy_df
  }
  ###### Sigflow/ Sigminer
  if (sigflow){
    fit_legacy <-  sig_fit(
      catalogue_matrix = t(mut_mat_input),
      sig_index = "ALL",
      sig_db = "legacy",
      mode = "SBS",
      method = "QP",
      return_class = "data.table",
      return_error = TRUE
    )
    output_fit(fit_legacy, result_dir = "sigflow" , mut_type = "legacy")
    sigflow_legacy_df = as_tibble(read.csv("sigflow/legacy_fitting_reconstruction_errors.csv"))%>%
      mutate(toolname = "sigflow")
    if (!exists("final_legacy_df")){
        final_legacy_df = sigflow_legacy_df[FALSE,]
    }
    final_legacy_df = rbind(final_legacy_df ,sigflow_legacy_df)    
  }
  ########### DeconstructSigs
  if (deconstructSigs){
    v2_input_matrices_list <- mut_mat_input
    v2_input_matrices_list <- apply(v2_input_matrices_list, 1, as.data.frame)
    v2_input_matrices_list <- lapply(v2_input_matrices_list, t)
    v2_input_matrices_list <- lapply(v2_input_matrices_list, as.data.frame)
    v2_signatures_list <- lapply(v2_input_matrices_list,whichSignatures,contexts.needed = TRUE,signatures.ref = signatures.cosmic)
    v2_exposures_list <-lapply(v2_signatures_list, function(x) x$weights)
    v2_diff_list <- lapply(v2_signatures_list, function(x) x$diff)
    v2_error_list <- lapply(seq(length(v2_diff_list))  , function(i)FrobeniusNorm_deconstructsigs(i,v2_input_matrices_list,v2_exposures_list,signatures.cosmic))
    dir.create("./deconstructsigs_results")
    write.csv(unlist(v2_error_list),row.names = rownames(mut_mat_input),"./deconstructsigs_results/legacy_sample_errors.csv")
    write.csv(do.call(rbind, lapply(v2_exposures_list, data.frame)),row.names = rownames(mut_mat_input),"./deconstructsigs_results/legacy_sample_exposures.csv")
    decon_legacy_df = as_tibble(read.csv("deconstructsigs_results/legacy_sample_errors.csv")) %>% 
      rename( sample = X,  error = x) %>%
      mutate(toolname = "deconstructsigs")
    if (!exists("final_legacy_df")){
        final_legacy_df = decon_legacy_df[FALSE,]
    }
    final_legacy_df = rbind(final_legacy_df , decon_legacy_df)    
  }
  ######### Sigfit
  if (sigfit){
    mut_mat = mut_mat_input
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
    errors_per_sample <- sapply(seq(ncol(expo_mat)), function(i) FrobeniusNorm(catalogue_matrix[i,], sig_matrix, expo_mat[, i]))
    rmse <- sqrt(mean(errors_per_sample ^ 2))
    errors_per_sample <- as.data.frame(errors_per_sample)
    row.names(errors_per_sample) <- row.names(sample_exposures$mean)
    dir.create("./sigfit_results")
    write.csv(errors_per_sample,"./sigfit_results/legacy_sample_errors.csv")
    write.csv(sample_exposures$mean, "./sigfit_results/legacy_sample_exposures.csv")
    sigfit_legacy_df = as_tibble(read.csv("sigfit_results/legacy_sample_errors.csv"))%>% 
      rename( sample = X,  error = errors_per_sample)%>%
      mutate(toolname = "sigfit")
    if (!exists("final_legacy_df")){
        final_legacy_df = sigfit_legacy_df[FALSE,]
    }
    final_legacy_df = rbind(final_legacy_df , sigfit_legacy_df)    
  }
  final_legacy_df = final_legacy_df %>% 
    mutate(mode = "Legacy SBS")
  return(final_legacy_df)
}
#########################################
# SBS
#########################################
runSBS = function(ref_genome = "hg19",
                  mut_mat_input, 
                  sbs_dropped_signatures = c("SBS27","SBS43","SBS45","SBS46","SBS47","SBS48","SBS49","SBS50","SBS51","SBS52","SBS53","SBS54","SBS55","SBS56","SBS57","SBS58","SBS59","SBS60"),
                  mutationalPatterns = TRUE, sigflow = TRUE , sigfit = TRUE , deconstructSigs = TRUE) {
  mut_mat = mut_mat_input

  ####### MutationalPatterns 
  if (mutationalPatterns){
    mut_mat <- t(mut_mat)
    ref_genome <- ref_genome
    # COSMIC 2019 SBS signatures
    # Strict signature fitting iteratively reduces the number of signatures used for refitting by removing the signature with least contribution
    signatures = t(SBS_signatures)
    signatures <- signatures[!row.names(signatures)%in%sbs_dropped_signatures,]
    signatures = as.matrix(t(signatures))
    rownames(signatures) = NULL
    strict_refit <- fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.004)
    fit_res_strict <- strict_refit$fit_res
    mut_patterns_strict_errors <-    as.matrix.data.frame(t(mut_mat - fit_res_strict$reconstructed))
    mut_patterns_strict_errors <-    apply(mut_patterns_strict_errors  ,1,FrobeniusNorm_mutational_patterns)
    write.csv(t(fit_res_strict$contribution),"./mutational_patterns_results/strict_sample_exposures.csv")
    write.csv(mut_patterns_strict_errors,row.names = colnames(mut_mat),"./mutational_patterns_results/strict_sample_errors.csv")
    # This is traditional decomposition
    fit_res <- fit_to_signatures(mut_mat, signatures)
    mut_patterns_errors <- as.matrix.data.frame(t(mut_mat - fit_res$reconstructed))
    mut_patterns_errors <-    apply(mut_patterns_errors  , 1, FrobeniusNorm_mutational_patterns)
    write.csv(t(fit_res$contribution) , "./mutational_patterns_results/sample_exposures.csv")
    write.csv(mut_patterns_errors, row.names = colnames(mut_mat), "./mutational_patterns_results/sample_errors.csv")
    mut_sbs_df = as_tibble(read.csv("mutational_patterns_results/sample_errors.csv") )%>% 
      rename( sample = X,  error = x)%>%
      mutate(toolname = "mutationalPatterns")
    final_sbs_df = mut_sbs_df
  }
  ###### Sigflow/ Sigminer
  if (sigflow){
    sbs_sigflow_signatures = t(SBS_signatures)
    sbs_sigflow_signatures <- sbs_sigflow_signatures[!row.names(sbs_sigflow_signatures)%in%sbs_dropped_signatures,]
    sbs_sigflow_signatures = t(sbs_sigflow_signatures)
    
    fit_sbs <-  sig_fit(
      catalogue_matrix = t(mut_mat_input),
      sig= sbs_sigflow_signatures, 
      sig_db = "SBS",
      mode = "SBS",
      method = "QP",
      return_class = "data.table",
      return_error = TRUE
    )
    output_fit(fit_sbs, result_dir =  "sigflow" , mut_type = "SBS")
    sigflow_sbs_df = as_tibble(read.csv("sigflow/SBS_fitting_reconstruction_errors.csv"))%>%
      mutate(toolname = "sigflow")
    if (!exists("final_sbs_df")){
        final_sbs_df = sigflow_sbs_df[FALSE,]
    }
    final_sbs_df = rbind(final_sbs_df, sigflow_sbs_df)
  }  
  ########### DeconstructSigs
  if (deconstructSigs){
    sbs_signatures = t(SBS_signatures)
    sbs_signatures = as.data.frame(sbs_signatures)
    sbs_signatures <- sbs_signatures[!row.names(sbs_signatures)%in%sbs_dropped_signatures,]
    v2_input_matrices_list <- mut_mat_input
    v2_input_matrices_list <- apply(v2_input_matrices_list, 1, as.data.frame)
    v2_input_matrices_list <- lapply(v2_input_matrices_list, t)
    v2_input_matrices_list <- lapply(v2_input_matrices_list, as.data.frame)
    v3_input_matrices_list <- v2_input_matrices_list
    v3_signatures_list <-  lapply(v3_input_matrices_list,whichSignatures,contexts.needed = TRUE,signatures.ref = sbs_signatures)
    v3_exposures_list <- lapply(v3_signatures_list, function(x)x$weights)
    v3_diff_list <- lapply(v3_signatures_list, function(x) x$diff)
    v3_error_list <- lapply(seq(length(v3_diff_list))  , function(i) FrobeniusNorm_deconstructsigs(i,v3_input_matrices_list,v3_exposures_list,sbs_signatures))
    write.csv(unlist(v3_error_list),row.names = rownames(mut_mat_input),"./deconstructsigs_results/sbs_sample_errors.csv")
    write.csv(do.call(rbind, lapply(v3_exposures_list, data.frame)),row.names = rownames(mut_mat_input),"./deconstructsigs_results/sbs_sample_exposures.csv")
    decon_sbs_df = as_tibble(read.csv("deconstructsigs_results/sbs_sample_errors.csv")) %>% 
      rename( sample = X,  error = x) %>%
      mutate(toolname = "deconstructsigs")
    if (!exists("final_sbs_df")){
        final_sbs_df = decon_sbs_df[FALSE,]
    }

    final_sbs_df = rbind(final_sbs_df, decon_sbs_df)
  }

  ######### Sigfit
  if (sigfit){
    sbs_sigfit_signatures = t(SBS_signatures)
    sbs_sigfit_signatures <- sbs_sigfit_signatures[!row.names(sbs_sigfit_signatures)%in%sbs_dropped_signatures,]
    samples_fit_sbs <- fit_signatures(counts = mut_mat_input ,signatures = sbs_sigfit_signatures,iter = 3000,warmup = 1000,chains = 1,seed = 1756)
    sample_exposures_sbs <-   retrieve_pars(samples_fit_sbs, "exposures")
    catalogue_matrix_sbs <- t(as.data.frame(mut_mat))
    sig_matrix_sbs <-  t(as.matrix(samples_fit_sbs[["data"]][["signatures"]]))
    expo_mat_sbs <- t(as.matrix(sample_exposures_sbs$mean))
    errors_per_sample_sbs <-    sapply(seq(ncol(expo_mat_sbs)), function(i) FrobeniusNorm(catalogue_matrix_sbs[i,], sig_matrix_sbs, expo_mat_sbs[, i]))
    rmse_sbs <- sqrt(mean(errors_per_sample_sbs ^ 2))
    errors_per_sample_sbs <- as.data.frame(errors_per_sample_sbs)
    row.names(errors_per_sample_sbs) <- row.names(sample_exposures_sbs$mean)
    write.csv(errors_per_sample_sbs,"./sigfit_results/sbs_sample_errors.csv")
    write.csv(sample_exposures_sbs$mean, "./sigfit_results/sbs_sample_exposures.csv")
    sigfit_sbs_df = as_tibble(read.csv("sigfit_results/sbs_sample_errors.csv"))%>% 
      rename( sample = X,  error = errors_per_sample_sbs)%>%
      mutate(toolname = "sigfit")
    if (!exists("final_sbs_df")){
        final_sbs_df = sigfit_sbs_df[FALSE,]
    }
    final_sbs_df = rbind(final_sbs_df, sigfit_sbs_df)
  }
  
  final_sbs_df = final_sbs_df %>%
    mutate(mode = "V3 SBS")
  return(final_sbs_df)
}
#########################################
# Indels
#########################################
runIndel = function(ref_genome = "hg19",
                      mut_mat_id_input, mutationalPatterns = TRUE, sigflow = TRUE , sigfit = TRUE , deconstructSigs = TRUE) {
  ####### MutationalPatterns 
  if(mutationalPatterns){
    mut_mat_id = t(mut_mat_id_input)
    signatures_id = get_known_signatures(source = "COSMIC", muttype = "indel")
    strict_refit <- fit_to_signatures(mut_mat_id, signatures_id)
    mut_patterns_id_errors <- as.matrix.data.frame(t(mut_mat_id - strict_refit$reconstructed))
    mut_patterns_id_errors <-apply(mut_patterns_id_errors  , 1,FrobeniusNorm_mutational_patterns)
    write.csv(t(strict_refit$contribution),"./mutational_patterns_results/id_sample_exposures.csv")
    write.csv(mut_patterns_id_errors,row.names = colnames(mut_mat_id),"./mutational_patterns_results/id_sample_errors.csv"  )
    mut_id_df = as_tibble(read.csv("mutational_patterns_results/id_sample_errors.csv") )%>% 
      rename( sample = X,  error = x)%>%
      mutate(toolname = "mutationalPatterns")
    final_id_df = mut_id_df
  }
  ###### Sigflow/ Sigminer
  if (sigflow){
    fit_indel <-  sig_fit(
      catalogue_matrix = t(mut_mat_id_input),
      sig_index = "ALL",
      sig_db = "ID",
      mode = "ID",
      method = "QP",
      return_class = "data.table",
      return_error = TRUE
    )
    output_fit(fit_indel, result_dir =  "sigflow" , mut_type = "ID")
    sigflow_id_df = as_tibble(read.csv("sigflow/ID_fitting_reconstruction_errors.csv"))%>%
      mutate(toolname = "sigflow")
    if (!exists("final_id_df")){
        final_id_df = sigflow_id_df[FALSE,]
    }
    final_id_df = rbind(final_id_df, sigflow_id_df)
  }
  ########### DeconstructSigs
  if (deconstructSigs){
    indel_signatures = get_known_signatures(source = "COSMIC", muttype = "indel")
    mut_mat_id = mut_mat_id_input
    indel_mutational_profile <- data.frame(mut_mat_id)
    indel_mutational_profile <- indel_mutational_profile %>% rename_at(colnames(indel_mutational_profile) , ~ colnames(as.data.frame(t(indel_signatures))))
    indel_input_matrices_list <- lapply(1:dim(indel_mutational_profile)[1], function(x) as.data.frame(indel_mutational_profile[x[1], ]))
    indel_signatures_list <- lapply(indel_input_matrices_list,whichSignaturesLocal,contexts.needed = TRUE,signatures.ref = as.data.frame(t(indel_signatures)) )
    indel_exposures_list <-lapply(indel_signatures_list, function(x) x$weights)
    indel_diff_list <-lapply(indel_signatures_list, function(x) x$diff)
    indel_error_list <-lapply(seq(length(indel_diff_list))  , function(i) FrobeniusNorm_deconstructsigs(i,indel_input_matrices_list,indel_exposures_list,as.data.frame(t(indel_signatures))))
    write.csv(unlist(indel_error_list),row.names = rownames(mut_mat_id),"./deconstructsigs_results/indel_sample_errors.csv")
    write.csv(do.call(rbind, lapply(indel_exposures_list, data.frame)),"./deconstructsigs_results/indel_sample_exposures.csv")
    decon_id_df = as_tibble(read.csv("deconstructsigs_results/indel_sample_errors.csv")) %>% 
      rename( sample = X,  error = x) %>%
      mutate(toolname = "deconstructsigs")
    if (!exists("final_id_df")){
        final_id_df = decon_id_df[FALSE,]
    }
    final_id_df = rbind(final_id_df, decon_id_df)
  }
  ######### Sigfit
  if (sigfit){
    mut_mat_id = mut_mat_id_input
    indel_mutational_profile = mut_mat_id
    indel_signatures = get_known_signatures(source = "COSMIC", muttype = "indel")
    indel_samples_fit <- fit_signatures(counts = indel_mutational_profile ,signatures = as.data.frame(t(ID_signatures)),iter = 3000,warmup = 1000,chains = 1,seed = 1756)
    indel_sample_exposures <-retrieve_pars(indel_samples_fit, "exposures")
    sig_matrix <-t(as.matrix(indel_samples_fit[["data"]][["signatures"]]))
    expo_mat <- t(as.matrix(indel_sample_exposures$mean))
    errors_per_sample <- sapply(seq(ncol(expo_mat)), function(i) FrobeniusNorm(indel_mutational_profile[i,], sig_matrix, expo_mat[, i]))
    errors_per_sample <- as.data.frame(errors_per_sample)
    row.names(errors_per_sample) <-row.names(indel_sample_exposures$mean)
    write.csv(errors_per_sample,"./sigfit_results/indel_sample_errors.csv")
    write.csv(indel_sample_exposures$mean,"./sigfit_results/indel_sample_exposures.csv")
    sigfit_id_df = as_tibble(read.csv("sigfit_results/indel_sample_errors.csv"))%>% 
      rename( sample = X,  error = errors_per_sample)%>%
      mutate(toolname = "sigfit")
    final_id_df = rbind(final_id_df, sigfit_id_df)
  }
  if (!exists("final_id_df")){
      final_id_df = sigfit_id_df[FALSE,]
  }
  final_id_df = final_id_df %>%
    mutate(mode = "ID")
  return(final_id_df)  
}
#########################################
# DBS
#########################################
runDBS = function(ref_genome = "hg19",
                      mut_mat_dbs_input, mutationalPatterns = TRUE, sigflow = TRUE , sigfit = TRUE , deconstructSigs = TRUE) {

  mut_mat_dbs = mut_mat_dbs_input
  ####### MutationalPatterns 
  if (mutationalPatterns){
    mut_mat_dbs = t(mut_mat_dbs)
    signatures_dbs = get_known_signatures(source = "COSMIC", muttype = "dbs")
    strict_refit <- fit_to_signatures(mut_mat_dbs, signatures_dbs)
    mut_patterns_dbs_errors <- as.matrix.data.frame(t(mut_mat_dbs - strict_refit$reconstructed))
    mut_patterns_dbs_errors <-apply(mut_patterns_dbs_errors, 1,FrobeniusNorm_mutational_patterns)
    write.csv(t(strict_refit$contribution),"./mutational_patterns_results/dbs_sample_exposures.csv")
    write.csv(mut_patterns_dbs_errors,row.names = colnames(mut_mat_dbs),"./mutational_patterns_results/dbs_sample_errors.csv")
    mut_dbs_df = as_tibble(read.csv("mutational_patterns_results/dbs_sample_errors.csv") )%>% 
      rename( sample = X,  error = x)%>%
      mutate(toolname = "mutationalPatterns")
    final_dbs_df = mut_dbs_df
  }
  ###### Sigflow/ Sigminer
  if (sigflow){
    fit_dbs <-  sig_fit(
      catalogue_matrix = t(mut_mat_dbs_input),
      sig_index = "ALL",
      sig_db = "DBS",
      mode = "DBS",
      method = "QP",
      return_class = "data.table",
      return_error = TRUE
    )
    output_fit(fit_dbs, result_dir =  "sigflow" , mut_type = "DBS")
    sigflow_dbs_df = as_tibble(read.csv("sigflow/DBS_fitting_reconstruction_errors.csv"))%>%
      mutate(toolname = "sigflow")
  if (!exists("final_dbs_df")){
      final_dbs_df = sigflow_dbs_df[FALSE,]
  }

    final_dbs_df = rbind(final_dbs_df , sigflow_dbs_df)
  }
  ########### DeconstructSigs
  if (deconstructSigs){
    dbs_signatures = get_known_signatures(source = "COSMIC", muttype = "dbs")
    mut_mat_dbs = mut_mat_dbs_input
    dbs_mutational_profile <- data.frame(mut_mat_dbs)
    dbs_mutational_profile <- dbs_mutational_profile %>% rename_at(colnames(dbs_mutational_profile) , ~ colnames(as.data.frame(t(dbs_signatures))))
    dbs_input_matrices_list <- lapply(1:dim(dbs_mutational_profile)[1], function(x) as.data.frame(dbs_mutational_profile[x[1], ]))
    dbs_signatures_list <-lapply(dbs_input_matrices_list,whichSignaturesLocal,contexts.needed = TRUE,signatures.ref = as.data.frame(t(dbs_signatures)))
    dbs_exposures_list <- lapply(dbs_signatures_list, function(x) x$weights)
    dbs_diff_list <- lapply(dbs_signatures_list, function(x) x$diff)
    dbs_error_list <- lapply(seq(length(dbs_diff_list))  , function(i) FrobeniusNorm_deconstructsigs(i,dbs_input_matrices_list,dbs_exposures_list,as.data.frame(t(dbs_signatures))))
    write.csv(unlist(dbs_error_list),row.names = rownames(mut_mat_dbs),"./deconstructsigs_results/dbs_sample_errors.csv")
    write.csv(do.call(rbind, lapply(dbs_exposures_list, data.frame)),"./deconstructsigs_results/dbs_sample_exposures.csv"  )
    decon_dbs_df = as_tibble(read.csv("deconstructsigs_results/dbs_sample_errors.csv")) %>% 
      rename( sample = X,  error = x) %>%
      mutate(toolname = "deconstructsigs")
  if (!exists("final_dbs_df")){
      final_dbs_df = decon_dbs_df[FALSE,]
  }

    final_dbs_df = rbind(final_dbs_df , decon_dbs_df)
  }
  ######### Sigfit
  if (sigfit){
    mut_mat_dbs = mut_mat_dbs_input
    dbs_mutational_profile = mut_mat_dbs
    dbs_signatures = get_known_signatures(source = "COSMIC", muttype = "dbs")
    dbs_samples_fit <- fit_signatures(counts = dbs_mutational_profile , signatures =  as.data.frame(t(dbs_signatures)),iter = 3000,warmup = 1000,chains = 1,seed = 1756)
    dbs_sample_exposures <- retrieve_pars(dbs_samples_fit, "exposures")
    sig_matrix <-t(as.matrix(dbs_samples_fit[["data"]][["signatures"]]))
    expo_mat <- t(as.matrix(dbs_sample_exposures$mean))
    errors_per_sample_dbs <- sapply(seq(ncol(expo_mat)), function(i) FrobeniusNorm(dbs_mutational_profile[i,], sig_matrix, expo_mat[, i]))
    errors_per_sample_dbs <- as.data.frame(errors_per_sample_dbs)
    row.names(errors_per_sample_dbs) <- row.names(dbs_sample_exposures$mean)
    write.csv(errors_per_sample_dbs,"./sigfit_results/dbs_sample_errors.csv")
    write.csv(dbs_sample_exposures$mean,"./sigfit_results/dbs_sample_exposures.csv")
    sigfit_dbs_df = as_tibble(read.csv("sigfit_results/dbs_sample_errors.csv"))%>% 
      rename( sample = X,  error = errors_per_sample_dbs)%>%
      mutate(toolname = "sigfit")
  if (!exists("final_dbs_df")){
      final_dbs_df = sigfit_dbs_df[FALSE,]
  }

    final_dbs_df = rbind(final_dbs_df , sigfit_dbs_df)
  }
  final_dbs_df = final_dbs_df %>%
    mutate(mode = "dbs")
}




#########################################
# Run commands
#########################################

# Parsing input
args = commandArgs(trailingOnly=TRUE)
input_directory <- args[1]
genome_build = args[2]

data("cosmic_signatures_v2")
if (genome_build == "GRCh37"){
SBS_signatures <- read.delim(sep = "," , colClasses=c("Type"="character"), "../data/COSMIC_v3.2_SBS_GRCh37.txt") %>% column_to_rownames(., var = "Type")
DBS_signatures <- read.delim(sep = "," , colClasses=c("Type"="character"), "../data/COSMIC_v3.2_DBS_GRCh37.txt") %>% column_to_rownames(., var = "Type")
}
if (genome_build == "GRCh38"){
  SBS_signatures <- read.delim(sep = "," , colClasses=c("Type"="character"), "../data/COSMIC_v3.2_SBS_GRCh38.txt") %>% column_to_rownames(., var = "Type")
  DBS_signatures <- read.delim(sep = "," , colClasses=c("Type"="character"),"../data/COSMIC_v3.2_DBS_GRCh38.txt") %>% column_to_rownames(., var = "Type")
}
if (genome_build == "mm9"){
  SBS_signatures <- read.delim(sep = "," , colClasses=c("Type"="character"), "../data/COSMIC_v3.2_SBS_mm9.txt") %>% column_to_rownames(., var = "Type")
  DBS_signatures <- read.delim(sep = "," , colClasses=c("Type"="character"),"../data/COSMIC_v3.2_DBS_mm9.txt") %>% column_to_rownames(., var = "Type")
}
if (genome_build == "mm10"){
  SBS_signatures <- read.delim(sep = "," , colClasses=c("Type"="character"), "../data/COSMIC_v3.2_SBS_mm10.txt") %>% column_to_rownames(., var = "Type")
  DBS_signatures <- read.delim(sep = "," , colClasses=c("Type"="character"), "../data/COSMIC_v3.2_DBS_mm10.txt") %>% column_to_rownames(., var = "Type")
}
if (genome_build == "rn6"){
  SBS_signatures <- read.delim(sep = "," , colClasses=c("Type"="character"), "../data/COSMIC_v3.2_SBS_rn6.txt") %>% column_to_rownames(., var = "Type")
  DBS_signatures <- read.delim(sep = "," , colClasses=c("Type"="character"), "../data/COSMIC_v3.2_DBS_rn6.txt") %>% column_to_rownames(., var = "Type")
}
ID_signatures <- read.delim(sep = "," , colClasses=c("Type"="character"), "../data/COSMIC_v3.2_ID_GRCh37.txt") %>% column_to_rownames(., var = "Type")
legacy_signatures <- readRDS("../data/legacy_signatures.RDs")$db



# genome_build = "GRCh37"

setwd(input_directory)

if ( file.exists("output/SBS/MetaMutationalSigs.SBS96.all")) {
  data_matrix_stratton =  as.data.frame(read.csv("output/SBS/MetaMutationalSigs.SBS96.all" ,sep = "\t"))
  result <- data_matrix_stratton[-1]
  row.names(result) <- data_matrix_stratton$MutationType 
  data_matrix_stratton = result
  run_sbs = TRUE
}else{
  run_sbs = FALSE
}

if ( file.exists("output/ID/MetaMutationalSigs.ID83.all")) {
  data_matrix_id_stratton =  as.data.frame(read.csv("output/ID/MetaMutationalSigs.ID83.all" ,sep = "\t"))
  result <- data_matrix_id_stratton[-1]
  row.names(result) <- data_matrix_id_stratton$MutationType 
  data_matrix_id_stratton = result
  run_indel = TRUE
}else {
  run_indel = FALSE
}

if ( file.exists("output/DBS/MetaMutationalSigs.DBS78.all")) {
  data_matrix_dbs_stratton =  as.data.frame(read.csv("output/DBS/MetaMutationalSigs.DBS78.all" ,sep = "\t"))
  result <- data_matrix_dbs_stratton[-1]
  row.names(result) <- data_matrix_dbs_stratton$MutationType 
  data_matrix_dbs_stratton = result
  run_dbs = TRUE
}else {
  run_dbs = FALSE
}

dir.create("MetaMutationalResults")
setwd("MetaMutationalResults")

input_mutationalPatterns = as.logical(args[3])
input_sigflow = as.logical(args[4]) 
input_sigfit = as.logical(args[5]) 
input_deconstructSigs = as.logical(args[6]) 

# input_mutationalPatterns = TRUE
# input_sigflow = FALSE 
# input_sigfit = FALSE 
# input_deconstructSigs = FALSE

if (run_sbs){
  final_legacy_df = runLegacy( mut_mat_input = t(data_matrix_stratton), mutationalPatterns = input_mutationalPatterns, sigflow = input_sigflow , sigfit = input_sigfit , deconstructSigs = input_deconstructSigs )
  final_sbs_df = runSBS( mut_mat_input = t(data_matrix_stratton) ,mutationalPatterns = input_mutationalPatterns, sigflow = input_sigflow , sigfit = input_sigfit , deconstructSigs = input_deconstructSigs)
  final_df =  rbind(final_legacy_df, final_sbs_df)
}
if (run_indel){
  final_id_df = runIndel( mut_mat_id_input = t(data_matrix_id_stratton) ,mutationalPatterns = input_mutationalPatterns, sigflow = input_sigflow , sigfit = input_sigfit , deconstructSigs = input_deconstructSigs)
    if (!exists("final_df")){
        final_df = final_id_df[FALSE,]
    }
  final_df =  rbind(final_df, final_id_df)

}
if (run_dbs){
  final_dbs_df = runDBS( mut_mat_dbs_input = t(data_matrix_dbs_stratton) , mutationalPatterns = input_mutationalPatterns, sigflow = input_sigflow , sigfit = input_sigfit , deconstructSigs = input_deconstructSigs)
    if (!exists("final_df")){
        final_df = final_dbs_df[FALSE,]
    }
  final_df =  rbind(final_df, final_dbs_df)

}

ggboxplot(final_df, x = "toolname",
          y = "error",
          combine = TRUE,
          color = "mode", palette = "jco",
          ylab = "error", 
          add = "jitter",
          add.params = list(size = 0.1, jitter = 0.2)  
)
ggsave("rmse_box_plot.svg")
