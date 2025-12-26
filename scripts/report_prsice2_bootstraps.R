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

compute_roc_curve <- function(prs, y01) {
  ok <- is.finite(prs) & !is.na(y01)
  prs <- prs[ok]
  y01 <- y01[ok]
  n1 <- sum(y01 == 1L)
  n0 <- sum(y01 == 0L)
  if (n1 < 1 || n0 < 1) return(NULL)

  ord <- order(prs, decreasing = TRUE)
  y <- y01[ord]
  tp <- cumsum(y == 1L)
  fp <- cumsum(y == 0L)

  tpr <- tp / n1
  fpr <- fp / n0

  # anchors
  data.frame(FPR = c(0, fpr, 1), TPR = c(0, tpr, 1))
}

open_png_headless_safe <- function(path, width_px, height_px, res_dpi) {
  if (.Platform$OS.type == "windows") {
    grDevices::png(path, width = width_px, height = height_px, res = res_dpi)
    return(TRUE)
  }
  if (capabilities("cairo")) {
    grDevices::png(path, width = width_px, height = height_px, res = res_dpi, type = "cairo")
    return(TRUE)
  }
  FALSE
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

read_prsice_summary_table <- function(path) {
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(
    read.delim(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) < 1) return(NULL)
  df
}

maybe_write_excel <- function(path, sheets) {
  # sheets: named list of data.frames
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    message("WARN: R package 'openxlsx' not available; skipping Excel: ", path)
    return(FALSE)
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  wb <- openxlsx::createWorkbook()
  for (nm in names(sheets)) {
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeDataTable(wb, nm, sheets[[nm]])
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  TRUE
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

  # Headless-safe PNG: prefer cairo if available (avoids X11 requirement)
  wrote_png <- FALSE
  wrote_png <- tryCatch({
    ok <- open_png_headless_safe(out_png, width_px = 1400, height_px = 600, res_dpi = 150)
    if (ok) {
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
      wrote_png <<- TRUE
    }
    wrote_png
  }, error = function(e) FALSE)
  if (!wrote_png) {
    message("WARN: Could not create PNG (likely headless without cairo): ", out_png)
  }

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
    mean_prs_control = numeric(0),
    mean_prs_case = numeric(0),
    delta_mean_prs = numeric(0),
    median_prs_control = numeric(0),
    median_prs_case = numeric(0),
    delta_median_prs = numeric(0),
    auc = numeric(0),
    prs_r2 = numeric(0),
    p = numeric(0),
    empirical_p = numeric(0),
    num_snp = integer(0),
    threshold = character(0),
    stringsAsFactors = FALSE
  )

  roc_curves <- list()
  summary_rows <- list()

  ok_reps <- 0
  skip_reps <- 0

  for (rep in seq_len(args$bootstraps)) {
    boot_dir <- file.path(args$prsice2_dir, paste0("bootstrap_", rep))
    best_path <- file.path(boot_dir, "prsice2.best")
    sum_path <- file.path(boot_dir, "prsice2.summary")
    fam_path <- file.path(args$split_dir, paste0("bootstrap_", rep), "test.fam")

    if (!file.exists(best_path)) {
      message("WARN: missing .best for bootstrap_", rep, ": ", best_path)
      skip_reps <- skip_reps + 1
      next
    }
    if (!file.exists(fam_path)) {
      message("WARN: missing test.fam for bootstrap_", rep, ": ", fam_path)
      skip_reps <- skip_reps + 1
      next
    }

    best <- tryCatch(read_best(best_path), error = function(e) {
      message("WARN: failed to read .best for bootstrap_", rep, ": ", conditionMessage(e))
      NULL
    })
    fam <- tryCatch(read_fam(fam_path), error = function(e) {
      message("WARN: failed to read test.fam for bootstrap_", rep, ": ", conditionMessage(e))
      NULL
    })
    if (is.null(best) || is.null(fam)) {
      skip_reps <- skip_reps + 1
      next
    }

    merged <- merge(best[, c("FID", "IID", "PRS")], fam[, c("FID", "IID", "PHENO")], by = c("FID", "IID"))
    if (nrow(merged) == 0) {
      message("WARN: no FID/IID overlap for bootstrap_", rep)
      skip_reps <- skip_reps + 1
      next
    }

    mapping <- infer_case_control(merged$PHENO, args$case_value, args$control_value)

    group <- rep(NA_character_, nrow(merged))
    group[merged$PHENO == mapping$control] <- "Control"
    group[merged$PHENO == mapping$case] <- "Case"

    keep <- !is.na(group) & is.finite(merged$PRS)
    merged <- merged[keep, , drop = FALSE]
    group <- group[keep]

    if (nrow(merged) == 0) {
      message("WARN: no usable case/control rows for bootstrap_", rep)
      skip_reps <- skip_reps + 1
      next
    }

    prs_case <- merged$PRS[group == "Case"]
    prs_control <- merged$PRS[group == "Control"]

    mean_control <- mean(prs_control, na.rm = TRUE)
    mean_case <- mean(prs_case, na.rm = TRUE)
    median_control <- stats::median(prs_control, na.rm = TRUE)
    median_case <- stats::median(prs_case, na.rm = TRUE)

    y01 <- ifelse(group == "Case", 1L, 0L)
    auc <- compute_auc(merged$PRS, y01)

    if (!file.exists(sum_path)) {
      message("WARN: missing .summary for bootstrap_", rep, ": ", sum_path)
    }
    summ <- read_prsice_summary(sum_path)

    sum_tab <- read_prsice_summary_table(sum_path)
    if (!is.null(sum_tab)) {
      sum_tab$bootstrap <- rep
      summary_rows[[length(summary_rows) + 1]] <- sum_tab
    }
    prs_r2 <- if (!is.null(summ)) summ$prs_r2 else NA_real_
    p <- if (!is.null(summ)) summ$p else NA_real_
    emp_p <- if (!is.null(summ)) summ$empirical_p else NA_real_
    num_snp <- if (!is.null(summ)) summ$num_snp else NA_integer_
    thr <- if (!is.null(summ)) summ$threshold else NA_character_

    # Per-bootstrap plots: create in bootstrap folder AND copy into report folder
    title <- paste0(args$target, " bootstrap_", rep)
    boot_png <- file.path(boot_dir, "prs_case_control.png")
    boot_pdf <- file.path(boot_dir, "prs_case_control.pdf")
    plot_case_control(prs_case, prs_control, title, boot_png, boot_pdf)

    rep_dir <- file.path(args$out_dir, "bootstraps")
    dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
    rep_png <- file.path(rep_dir, paste0("bootstrap_", rep, "_prs_case_control.png"))
    rep_pdf <- file.path(rep_dir, paste0("bootstrap_", rep, "_prs_case_control.pdf"))

    if (file.exists(boot_pdf)) file.copy(boot_pdf, rep_pdf, overwrite = TRUE)
    if (file.exists(boot_png)) file.copy(boot_png, rep_png, overwrite = TRUE)

    metrics <- rbind(
      metrics,
      data.frame(
        rep = rep,
        n_total = nrow(merged),
        n_control = sum(group == "Control"),
        n_case = sum(group == "Case"),
        mean_prs_control = mean_control,
        mean_prs_case = mean_case,
        delta_mean_prs = mean_case - mean_control,
        median_prs_control = median_control,
        median_prs_case = median_case,
        delta_median_prs = median_case - median_control,
        auc = auc,
        prs_r2 = prs_r2,
        p = p,
        empirical_p = emp_p,
        num_snp = num_snp,
        threshold = thr,
        stringsAsFactors = FALSE
      )
    )

    roc <- compute_roc_curve(merged$PRS, y01)
    if (!is.null(roc)) {
      roc$rep <- rep
      roc_curves[[length(roc_curves) + 1]] <- roc
    }

    ok_reps <- ok_reps + 1
  }

  if (nrow(metrics) == 0) {
    stop("No bootstraps could be summarized (missing .best/.fam or no usable phenotypes).")
  }

  # Write summary table
  out_tsv <- file.path(args$out_dir, "prsice2_bootstrap_summary.tsv")
  write.table(metrics, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

  # Merge PRSice2 .summary tables across bootstraps
  summary_all <- NULL
  if (length(summary_rows) > 0) {
    summary_all <- tryCatch(do.call(rbind, summary_rows), error = function(e) NULL)
  }
  summary_tsv <- file.path(args$out_dir, "prsice2_summary_all_bootstraps.tsv")
  if (!is.null(summary_all)) {
    write.table(summary_all, file = summary_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    message("WARN: No .summary tables could be merged; not writing: ", summary_tsv)
  }

  # Write Excel workbook if possible
  xlsx_path <- file.path(args$out_dir, "prsice2_report.xlsx")
  sheets <- list(metrics = metrics)
  if (!is.null(summary_all)) sheets$prsice2_summary_all <- summary_all
  maybe_write_excel(xlsx_path, sheets)

  # Aggregate plots
  # AUC across bootstraps (different colors)
  auc_png <- file.path(args$out_dir, "auc_across_bootstraps.png")
  auc_pdf <- file.path(args$out_dir, "auc_across_bootstraps.pdf")

  wrote_auc_png <- FALSE
  wrote_auc_png <- tryCatch({
    ok <- open_png_headless_safe(auc_png, width_px = 1400, height_px = 700, res_dpi = 150)
    if (!ok) return(FALSE)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); dev.off() }, add = TRUE)
    par(mar = c(7, 4, 4, 1) + 0.1)
    cols <- grDevices::rainbow(nrow(metrics))
    barplot(metrics$auc, names.arg = metrics$rep, col = cols, las = 2, ylim = c(0, 1),
            main = paste0(args$target, " AUC across bootstraps"), xlab = "Bootstrap", ylab = "AUC")
    TRUE
  }, error = function(e) FALSE)
  if (!wrote_auc_png) message("WARN: Could not create PNG (likely headless without cairo): ", auc_png)

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

  wrote_r2_png <- FALSE
  wrote_r2_png <- tryCatch({
    ok <- open_png_headless_safe(r2_png, width_px = 1400, height_px = 700, res_dpi = 150)
    if (!ok) return(FALSE)
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); dev.off() }, add = TRUE)
    par(mar = c(7, 4, 4, 1) + 0.1)
    cols <- grDevices::rainbow(nrow(metrics))
    barplot(metrics$prs_r2, names.arg = metrics$rep, col = cols, las = 2,
            main = paste0(args$target, " PRS.R2 across bootstraps"), xlab = "Bootstrap", ylab = "PRS.R2")
    TRUE
  }, error = function(e) FALSE)
  if (!wrote_r2_png) message("WARN: Could not create PNG (likely headless without cairo): ", r2_png)

  # Case vs control PRS (means) across bootstraps
  mean_png <- file.path(args$out_dir, "mean_prs_case_control_across_bootstraps.png")
  mean_pdf <- file.path(args$out_dir, "mean_prs_case_control_across_bootstraps.pdf")

  plot_mean_case_control <- function(open_device) {
    open_device()
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); dev.off() }, add = TRUE)
    par(mar = c(7, 4, 4, 1) + 0.1)

    x <- seq_along(metrics$rep)
    cols <- grDevices::rainbow(length(x))
    y_all <- range(c(metrics$mean_prs_control, metrics$mean_prs_case), finite = TRUE)
    plot(x, metrics$mean_prs_control, type = "b", pch = 16, col = "darkgreen", xaxt = "n",
         xlab = "Bootstrap", ylab = "Mean PRS", ylim = y_all,
         main = paste0(args$target, " mean PRS: Control vs Case"))
    axis(1, at = x, labels = metrics$rep, las = 2)
    lines(x, metrics$mean_prs_case, type = "b", pch = 16, col = "darkred")
    legend("topleft", legend = c("Control", "Case"), col = c("darkgreen", "darkred"), lwd = 2, pch = 16, bty = "n")
  }

  wrote_mean_png <- tryCatch({
    ok <- open_png_headless_safe(mean_png, width_px = 1400, height_px = 700, res_dpi = 150)
    if (!ok) return(FALSE)
    plot_mean_case_control(function() {})
    TRUE
  }, error = function(e) FALSE)
  if (!wrote_mean_png) message("WARN: Could not create PNG (likely headless without cairo): ", mean_png)
  plot_mean_case_control(function() grDevices::pdf(mean_pdf, width = 12, height = 6, useDingbats = FALSE))

  # Delta mean PRS across bootstraps
  d_png <- file.path(args$out_dir, "delta_mean_prs_case_minus_control.png")
  d_pdf <- file.path(args$out_dir, "delta_mean_prs_case_minus_control.pdf")

  plot_delta <- function(open_device) {
    open_device()
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); dev.off() }, add = TRUE)
    par(mar = c(7, 4, 4, 1) + 0.1)
    cols <- grDevices::rainbow(nrow(metrics))
    barplot(metrics$delta_mean_prs, names.arg = metrics$rep, col = cols, las = 2,
            main = paste0(args$target, " mean(PRS_case) - mean(PRS_control)"),
            xlab = "Bootstrap", ylab = "Delta mean PRS")
    abline(h = 0, lty = 2)
  }

  wrote_d_png <- tryCatch({
    ok <- open_png_headless_safe(d_png, width_px = 1400, height_px = 700, res_dpi = 150)
    if (!ok) return(FALSE)
    plot_delta(function() {})
    TRUE
  }, error = function(e) FALSE)
  if (!wrote_d_png) message("WARN: Could not create PNG (likely headless without cairo): ", d_png)
  plot_delta(function() grDevices::pdf(d_pdf, width = 12, height = 6, useDingbats = FALSE))

  # ROC plot across bootstraps (all curves on one plot)
  roc_png <- file.path(args$out_dir, "roc_across_bootstraps.png")
  roc_pdf <- file.path(args$out_dir, "roc_across_bootstraps.pdf")
  roc_df <- if (length(roc_curves) > 0) do.call(rbind, roc_curves) else NULL

  plot_roc <- function(open_device) {
    open_device()
    op <- par(no.readonly = TRUE)
    on.exit({ par(op); dev.off() }, add = TRUE)
    par(mar = c(5, 5, 4, 8) + 0.1, xpd = TRUE)

    plot(c(0, 1), c(0, 1), type = "n", xlab = "False Positive Rate", ylab = "True Positive Rate",
         main = paste0(args$target, " ROC across bootstraps"), xaxs = "i", yaxs = "i")
    abline(0, 1, lty = 2, col = "black")

    if (!is.null(roc_df)) {
      reps <- sort(unique(roc_df$rep))
      cols <- grDevices::rainbow(length(reps))
      for (i in seq_along(reps)) {
        r <- reps[[i]]
        sub <- roc_df[roc_df$rep == r, , drop = FALSE]
        lines(sub$FPR, sub$TPR, col = cols[[i]], lwd = 2)
      }
      auc_map <- metrics$auc
      names(auc_map) <- metrics$rep
      leg <- paste0("bootstrap_", reps, " (AUC=", sprintf("%.3f", auc_map[as.character(reps)]), ")")
      legend("right", inset = c(-0.25, 0), legend = leg, col = cols, lwd = 2, cex = 0.75, bty = "n")
    } else {
      text(0.5, 0.5, "No ROC curves available", cex = 1.2)
    }
  }

  wrote_roc_png <- tryCatch({
    ok <- open_png_headless_safe(roc_png, width_px = 1400, height_px = 900, res_dpi = 150)
    if (!ok) return(FALSE)
    plot_roc(function() {})
    TRUE
  }, error = function(e) FALSE)
  if (!wrote_roc_png) message("WARN: Could not create PNG (likely headless without cairo): ", roc_png)

  plot_roc(function() grDevices::pdf(roc_pdf, width = 12, height = 8, useDingbats = FALSE))

  grDevices::pdf(r2_pdf, width = 12, height = 6, useDingbats = FALSE)
  op <- par(no.readonly = TRUE)
  par(mar = c(7, 4, 4, 1) + 0.1)
  cols <- grDevices::rainbow(nrow(metrics))
  barplot(metrics$prs_r2, names.arg = metrics$rep, col = cols, las = 2,
          main = paste0(args$target, " PRS.R2 across bootstraps"), xlab = "Bootstrap", ylab = "PRS.R2")
  par(op)
  dev.off()

  message("Bootstraps summarized: ", ok_reps, "; skipped: ", skip_reps)
  message("Wrote summary: ", out_tsv)
  if (!is.null(summary_all)) message("Wrote merged summaries: ", summary_tsv)
  message("Excel (if available): ", xlsx_path)
  message("Report folder: ", args$out_dir)
  message("Aggregate plots (PDF): ", auc_pdf, ", ", r2_pdf, ", ", roc_pdf)
  message("Extra comparison plots (PDF): ", mean_pdf, ", ", d_pdf)
}

main()
