#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- function(exit_code = 1) {
  cat(
    "Usage:\n",
    "  Rscript scripts/format_gwas_for_prsice.R --in GWAS.assoc.logistic --bim target.bim --out base.tsv\n\n",
    "Input (PLINK .assoc.logistic): whitespace-delimited, must contain columns: CHR SNP BP A1 TEST OR P\n",
    "Output (PRSice base): tab-delimited columns: SNP\tCHR\tBP\tA1\tA2\tOR\tP\n",
    sep = ""
  )
  quit(status = exit_code)
}

get_flag_value <- function(flag) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(NULL)
  if (idx[length(idx)] == length(args)) {
    cat(sprintf("ERROR: Missing value for %s\n", flag), file = stderr())
    usage(2)
  }
  args[idx[length(idx)] + 1]
}

in_path <- get_flag_value("--in")
bim_path <- get_flag_value("--bim")
out_path <- get_flag_value("--out")

test_filter <- get_flag_value("--test")
if (is.null(test_filter) || test_filter == "") test_filter <- "ADD"

if (is.null(in_path) || is.null(bim_path) || is.null(out_path)) {
  cat("ERROR: --in, --bim, and --out are required\n", file = stderr())
  usage(2)
}

if (!file.exists(in_path)) {
  cat(sprintf("ERROR: Input file not found: %s\n", in_path), file = stderr())
  quit(status = 2)
}
if (!file.exists(bim_path)) {
  cat(sprintf("ERROR: BIM file not found: %s\n", bim_path), file = stderr())
  quit(status = 2)
}

have_dt <- requireNamespace("data.table", quietly = TRUE)

cat(sprintf("Reading assoc: %s\n", in_path), file = stderr())
gwas <- NULL
if (have_dt) {
  gwas <- tryCatch(
    data.table::fread(in_path, header = TRUE, data.table = TRUE, showProgress = FALSE),
    error = function(e) NULL
  )
}
if (is.null(gwas)) {
  # Base R fallback (whitespace-delimited)
  gwas <- tryCatch(
    read.table(in_path, header = TRUE, sep = "", stringsAsFactors = FALSE, fill = TRUE, comment.char = ""),
    error = function(e) {
      cat(sprintf("ERROR: failed to read assoc with base R: %s\n", e$message), file = stderr())
      NULL
    }
  )
}

if (is.null(gwas) || ncol(gwas) == 0) {
  cat("ERROR: Failed to read assoc file\n", file = stderr())
  quit(status = 2)
}

cat(sprintf("Reading bim:   %s\n", bim_path), file = stderr())
bim <- NULL
if (have_dt) {
  bim <- tryCatch(
    data.table::fread(bim_path, header = FALSE, data.table = TRUE, showProgress = FALSE),
    error = function(e) NULL
  )
}
if (is.null(bim)) {
  bim <- tryCatch(
    read.table(bim_path, header = FALSE, sep = "", stringsAsFactors = FALSE, fill = TRUE, comment.char = ""),
    error = function(e) {
      cat(sprintf("ERROR: failed to read bim with base R: %s\n", e$message), file = stderr())
      NULL
    }
  )
}

if (is.null(bim) || ncol(bim) < 6) {
  cat("ERROR: BIM file must have at least 6 columns\n", file = stderr())
  quit(status = 2)
}

if (have_dt) {
  data.table::setDT(bim)
  data.table::setnames(bim, c("CHR_BIM", "SNP", "CM", "BP_BIM", "A1_BIM", "A2_BIM"))
} else {
  colnames(bim)[1:6] <- c("CHR_BIM", "SNP", "CM", "BP_BIM", "A1_BIM", "A2_BIM")
}

req <- c("CHR", "SNP", "BP", "A1", "TEST", "OR", "P")
missing <- req[!req %in% names(gwas)]
if (length(missing) > 0) {
  cat("ERROR: assoc header missing required columns:\n", file = stderr())
  cat(paste0("  - ", missing, collapse = "\n"), "\n", file = stderr())
  cat("Have:\n", file = stderr())
  cat(paste0("  ", names(gwas), collapse = "\n"), "\n", file = stderr())
  quit(status = 2)
}

# Filter to ADD test
if (!is.null(test_filter) && test_filter != "") {
  gwas <- gwas[gwas$TEST == test_filter, , drop = FALSE]
}

# Join BIM alleles and compute A2 as the non-effect allele.
if (have_dt) {
  data.table::setDT(gwas)
  data.table::setkey(bim, SNP)
  data.table::setkey(gwas, SNP)
  merged <- merge(gwas, bim[, .(SNP, A1_BIM, A2_BIM)], by = "SNP", all.x = TRUE)

  merged[, A1_UP := toupper(as.character(A1))]
  merged[, A1_BIM_UP := toupper(as.character(A1_BIM))]
  merged[, A2_BIM_UP := toupper(as.character(A2_BIM))]

  merged[, A2 := data.table::fifelse(
    !is.na(A1_UP) & A1_UP == A1_BIM_UP, A2_BIM_UP,
    data.table::fifelse(!is.na(A1_UP) & A1_UP == A2_BIM_UP, A1_BIM_UP, NA_character_)
  )]
} else {
  merged <- merge(gwas, bim[, c("SNP", "A1_BIM", "A2_BIM")], by = "SNP", all.x = TRUE)
  merged$A1_UP <- toupper(as.character(merged$A1))
  merged$A1_BIM_UP <- toupper(as.character(merged$A1_BIM))
  merged$A2_BIM_UP <- toupper(as.character(merged$A2_BIM))
  merged$A2 <- ifelse(
    !is.na(merged$A1_UP) & merged$A1_UP == merged$A1_BIM_UP,
    merged$A2_BIM_UP,
    ifelse(!is.na(merged$A1_UP) & merged$A1_UP == merged$A2_BIM_UP, merged$A1_BIM_UP, NA)
  )
}

if (have_dt) {
  out <- merged[, .(
    SNP = as.character(SNP),
    CHR = as.character(CHR),
    BP = as.integer(BP),
    A1 = toupper(as.character(A1)),
    A2 = toupper(as.character(A2)),
    OR = as.numeric(OR),
    P = as.numeric(P)
  )]
} else {
  out <- data.frame(
    SNP = as.character(merged$SNP),
    CHR = as.character(merged$CHR),
    BP = as.integer(merged$BP),
    A1 = toupper(as.character(merged$A1)),
    A2 = toupper(as.character(merged$A2)),
    OR = as.numeric(merged$OR),
    P = as.numeric(merged$P),
    stringsAsFactors = FALSE
  )
}

# Drop invalid rows
if (have_dt) {
  out <- out[!is.na(SNP) & SNP != "" & !is.na(CHR) & CHR != "" & !is.na(BP) & !is.na(A1) & A1 != "" & !is.na(OR) & !is.na(P)]
} else {
  keep <- !is.na(out$SNP) & out$SNP != "" & !is.na(out$CHR) & out$CHR != "" & !is.na(out$BP) & !is.na(out$A1) & out$A1 != "" & !is.na(out$OR) & !is.na(out$P)
  out <- out[keep, , drop = FALSE]
}

# Ensure output dir exists
out_dir <- dirname(out_path)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

cat(sprintf("Writing base: %s\n", out_path), file = stderr())
if (have_dt) {
  data.table::fwrite(out, file = out_path, sep = "\t", quote = FALSE, col.names = TRUE)
} else {
  write.table(out, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}
cat(sprintf("Rows written: %d\n", nrow(out)), file = stderr())
