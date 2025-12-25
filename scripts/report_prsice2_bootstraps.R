#!/usr/bin/env Rscript

# Report PRSice2 results across bootstraps:
# - For each bootstrap: read PRSice2 .best and test .fam, compute AUC, and plot case/control PRS distribution
# - For all bootstraps: parse PRS.R2 from PRSice2 .summary and plot PRS.R2 + AUC across bootstraps
#
# Uses base R only (no external packages) for maximum portability on clusters.

parse_args <- function(args) {
  out <- list(
    results_dir = NULL,
    target = NULL,
    bootstraps = NULL,
    prsice2_dir = NULL,
    split_dir = NULL,
    out_dir = NULL,
    case_value = NA_integer_,
    control_value = NA_integer_
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("-h", "--help")) {
      cat(
        paste(
          "Summarize PRSice2 outputs across bootstraps.",
          "\n\nUsage:",
          "  Rscript scripts/report_prsice2_bootstraps.R \\",
          "    --results-dir results --target TARGET --bootstraps N \\",
          "    --prsice2-dir results/TARGET/prs/PRSice2 \\",
          "    --split-dir results/TARGET/split_train_pca_gwas \\",
          "    --out-dir results/TARGET/prs/PRSice2/report \\",
          "    [--case-value 2 --control-value 1]",
          sep = "\n"
        )
      )
      quit(status = 0)
    }

    if (key == "--results-dir") { out$results_dir <- args[[i + 1]]; i <- i + 2; next }
    if (key == "--target") { out$target <- args[[i + 1]]; i <- i + 2; next }
    if (key == "--bootstraps") { out$bootstraps <- as.integer(args[[i + 1]]); i <- i + 2; next }
    if (key == "--prsice2-dir") { out$prsice2_dir <- args[[i + 1]]; i <- i + 2; next }
    if (key == "--split-dir") { out$split_dir <- args[[i + 1]]; i <- i + 2; next }
    if (key == "--out-dir") { out$out_dir <- args[[i + 1]]; i <- i + 2; next }
    if (key == "--case-value") { out$case_value <- as.integer(args[[i + 1]]); i <- i + 2; next }
    if (key == "--control-value") { out$control_value <- as.integer(args[[i + 1]]); i <- i + 2; next }

    stop(paste0("Unknown argument: ", key))
  }

  if (is.null(out$results_dir) || is.null(out$target) || is.null(out$bootstraps) ||
      is.null(out$prsice2_dir) || is.null(out$split_dir) || is.null(out$out_dir)) {
    stop("Missing required args. Use --help.")
  }
  if (!is.finite(out$bootstraps) || out$bootstraps < 1) stop("--bootstraps must be >= 1")

  out
}

read_best <- function(path) {
  if (!file.exists(path)) stop(paste0("Missing .best: ", path))
  best <- read.table(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("FID", "IID") %in% names(best))) stop(paste0(".best missing FID/IID: ", path))
  if (!("PRS" %in% names(best))) {
    if ("SCORE" %in% names(best)) {
      best$PRS <- best$SCORE
    } else {
      stop(paste0(".best missing PRS/SCORE column: ", path))
    }
  }
  best$FID <- as.character(best$FID)
  best$IID <- as.character(best$IID)
  best$PRS <- suppressWarnings(as.numeric(best$PRS))
  best
}

read_fam <- function(path) {
  if (!file.exists(path)) stop(paste0("Missing .fam: ", path))
  fam <- read.table(path, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(fam) < 6) stop(paste0("Unexpected .fam format (<6 cols): ", path))
  fam <- fam[, 1:6]
  names(fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")
  fam$FID <- as.character(fam$FID)
  fam$IID <- as.character(fam$IID)
  fam$PHENO <- suppressWarnings(as.integer(fam$PHENO))
  fam
}

infer_case_control <- function(pheno, case_value, control_value) {
  if (!is.na(case_value) && !is.na(control_value)) {
    return(list(case = case_value, control = control_value))
  }

  vals <- sort(unique(pheno[!is.na(pheno)]))

  # Common encodings
  if (all(vals %in% c(-9L, 0L, 1L, 2L))) return(list(case = 2L, control = 1L))
  if (all(vals %in% c(0L, 1L))) return(list(case = 1L, control = 0L))

  # Fallback: treat max positive as case, min positive as control
  vals_pos <- vals[vals >= 0]
  if (length(vals_pos) >= 2) return(list(case = max(vals_pos), control = min(vals_pos)))

  # Last resort
  list(case = 2L, control = 1L)
}

compute_auc <- function(prs, is_case01) {
  ok <- is.finite(prs) & !is.na(is_case01)
  prs <- prs[ok]
  is_case01 <- is_case01[ok]
  n1 <- sum(is_case01 == 1L)
  n0 <- sum(is_case01 == 0L)
  if (n1 < 1 || n0 < 1) return(NA_real_)
  r <- rank(prs, ties.method = "average")
  sum_r_case <- sum(r[is_case01 == 1L])
  (sum_r_case - n1 * (n1 + 1) / 2) / (n1 * n0)
}

read_prsice_summary <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(
    read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) < 1) return(NULL)

  # Prefer Set==Base if present
  if ("Set" %in% names(df)) {
    idx <- which(df$Set == "Base")
    if (length(idx) >= 1) df <- df[idx[1], , drop = FALSE]
  }

  numify <- function(x) suppressWarnings(as.numeric(x))
  out <- list(
    threshold = if ("Threshold" %in% names(df)) as.character(df[[1, "Threshold"]]) else NA_character_,
    prs_r2 = if ("PRS.R2" %in% names(df)) numify(df[[1, "PRS.R2"]]) else NA_real_,
    p = if ("P" %in% names(df)) numify(df[[1, "P"]]) else NA_real_,
    num_snp = if ("Num_SNP" %in% names(df)) suppressWarnings(as.integer(df[[1, "Num_SNP"]])) else NA_integer_,
    empirical_p = if ("Empirical-P" %in% names(df)) numify(df[[1, "Empirical-P"]]) else NA_real_
  )
  out
}

plot_case_control <- function(prs_case, prs_control, title, out_png, out_pdf) {
  dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

  # 2-panel plot: boxplot + density
  make_plot <- function(device_fun) {
    device_fun()
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); dev.off() }, add = TRUE)

    par(mfrow = c(1, 2), mar = c(5, 4, 4, 1) + 0.1)

    boxplot(
      list(Control = prs_control, Case = prs_case),
      col = c("lightgreen", "lightcoral"),
      main = title,
      ylab = "PRS"
    )

    d0 <- density(prs_control, na.rm = TRUE)
    d1 <- density(prs_case, na.rm = TRUE)
    ylim <- range(c(d0$y, d1$y), finite = TRUE)
    xlim <- range(c(d0$x, d1$x), finite = TRUE)

    plot(d0, col = "darkgreen", lwd = 2, main = "PRS density", xlab = "PRS", ylim = ylim, xlim = xlim)
    lines(d1, col = "darkred", lwd = 2)
    legend("topright", legend = c("Control", "Case"), col = c("darkgreen", "darkred"), lwd = 2, bty = "n")
  }

  make_plot(function() grDevices::png(out_png, width = 1400, height = 600, res = 150))
  make_plot(function() grDevices::pdf(out_pdf, width = 12, height = 5, useDingbats = FALSE))
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

  metrics <- data.frame(
    rep = integer(0),
    n_total = integer(0),
    n_control = integer(0),
    n_case = integer(0),
    auc = numeric(0),
    prs_r2 = numeric(0),
    p = numeric(0),
    empirical_p = numeric(0),
    num_snp = integer(0),
    threshold = character(0),
    stringsAsFactors = FALSE
  )

  for (rep in seq_len(args$bootstraps)) {
    boot_dir <- file.path(args$prsice2_dir, paste0("bootstrap_", rep))
    best_path <- file.path(boot_dir, "prsice2.best")
    sum_path <- file.path(boot_dir, "prsice2.summary")
    fam_path <- file.path(args$split_dir, paste0("bootstrap_", rep), "test.fam")

    best <- read_best(best_path)
    fam <- read_fam(fam_path)

    merged <- merge(best[, c("FID", "IID", "PRS")], fam[, c("FID", "IID", "PHENO")], by = c("FID", "IID"))
    if (nrow(merged) == 0) stop(paste0("No FID/IID overlap for rep ", rep))

    mapping <- infer_case_control(merged$PHENO, args$case_value, args$control_value)

    group <- rep(NA_character_, nrow(merged))
    group[merged$PHENO == mapping$control] <- "Control"
    group[merged$PHENO == mapping$case] <- "Case"

    keep <- !is.na(group) & is.finite(merged$PRS)
    merged <- merged[keep, , drop = FALSE]
    group <- group[keep]

    if (nrow(merged) == 0) stop(paste0("No usable rows (case/control) for rep ", rep))

    prs_case <- merged$PRS[group == "Case"]
    prs_control <- merged$PRS[group == "Control"]

    y01 <- ifelse(group == "Case", 1L, 0L)
    auc <- compute_auc(merged$PRS, y01)

    summ <- read_prsice_summary(sum_path)
    prs_r2 <- if (!is.null(summ)) summ$prs_r2 else NA_real_
    p <- if (!is.null(summ)) summ$p else NA_real_
    emp_p <- if (!is.null(summ)) summ$empirical_p else NA_real_
    num_snp <- if (!is.null(summ)) summ$num_snp else NA_integer_
    thr <- if (!is.null(summ)) summ$threshold else NA_character_

    # Per-bootstrap plots written into each bootstrap folder
    title <- paste0(args$target, " bootstrap_", rep)
    out_png <- file.path(boot_dir, "prs_case_control.png")
    out_pdf <- file.path(boot_dir, "prs_case_control.pdf")
    plot_case_control(prs_case, prs_control, title, out_png, out_pdf)

    metrics <- rbind(
      metrics,
      data.frame(
        rep = rep,
        n_total = nrow(merged),
        n_control = sum(group == "Control"),
        n_case = sum(group == "Case"),
        auc = auc,
        prs_r2 = prs_r2,
        p = p,
        empirical_p = emp_p,
        num_snp = num_snp,
        threshold = thr,
        stringsAsFactors = FALSE
      )
    )
  }

  # Write summary table
  out_tsv <- file.path(args$out_dir, "prsice2_bootstrap_summary.tsv")
  write.table(metrics, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

  # Aggregate plots
  # AUC across bootstraps (different colors)
  auc_png <- file.path(args$out_dir, "auc_across_bootstraps.png")
  auc_pdf <- file.path(args$out_dir, "auc_across_bootstraps.pdf")

  grDevices::png(auc_png, width = 1400, height = 700, res = 150)
  op <- par(no.readonly = TRUE)
  par(mar = c(7, 4, 4, 1) + 0.1)
  cols <- grDevices::rainbow(nrow(metrics))
  barplot(metrics$auc, names.arg = metrics$rep, col = cols, las = 2, ylim = c(0, 1),
          main = paste0(args$target, " AUC across bootstraps"), xlab = "Bootstrap", ylab = "AUC")
  par(op)
  dev.off()

  grDevices::pdf(auc_pdf, width = 12, height = 6, useDingbats = FALSE)
  op <- par(no.readonly = TRUE)
  par(mar = c(7, 4, 4, 1) + 0.1)
  cols <- grDevices::rainbow(nrow(metrics))
  barplot(metrics$auc, names.arg = metrics$rep, col = cols, las = 2, ylim = c(0, 1),
          main = paste0(args$target, " AUC across bootstraps"), xlab = "Bootstrap", ylab = "AUC")
  par(op)
  dev.off()

  # PRS.R2 across bootstraps
  r2_png <- file.path(args$out_dir, "prs_r2_across_bootstraps.png")
  r2_pdf <- file.path(args$out_dir, "prs_r2_across_bootstraps.pdf")

  grDevices::png(r2_png, width = 1400, height = 700, res = 150)
  op <- par(no.readonly = TRUE)
  par(mar = c(7, 4, 4, 1) + 0.1)
  cols <- grDevices::rainbow(nrow(metrics))
  barplot(metrics$prs_r2, names.arg = metrics$rep, col = cols, las = 2,
          main = paste0(args$target, " PRS.R2 across bootstraps"), xlab = "Bootstrap", ylab = "PRS.R2")
  par(op)
  dev.off()

  grDevices::pdf(r2_pdf, width = 12, height = 6, useDingbats = FALSE)
  op <- par(no.readonly = TRUE)
  par(mar = c(7, 4, 4, 1) + 0.1)
  cols <- grDevices::rainbow(nrow(metrics))
  barplot(metrics$prs_r2, names.arg = metrics$rep, col = cols, las = 2,
          main = paste0(args$target, " PRS.R2 across bootstraps"), xlab = "Bootstrap", ylab = "PRS.R2")
  par(op)
  dev.off()

  message("Wrote summary: ", out_tsv)
  message("Wrote plots: ", auc_png, ", ", r2_png)
}

main()
