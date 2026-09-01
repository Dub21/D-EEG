#!/usr/bin/env Rscript
# Export des parametres normatifs (mu, sigma, nu, tau) sur une grille age x sexe,
# a partir des objets gamlss ajustes. A EXECUTER SUR LE CLUSTER.
#
# Sortie : une table compacte, agregee, sans aucune donnee individuelle.
# C'est cette table (et elle seule) qui part vers l'app / shinyapps.io.
#
# Usage :
#   Rscript export_normative_params.R <MODEL_DIR> <OUT_DIR> [N_AGE]

suppressPackageStartupMessages({ library(gamlss) })

args     <- commandArgs(trailingOnly = TRUE)
MODEL_DIR <- args[1]
OUT_DIR   <- args[2]
N_AGE     <- if (length(args) >= 3) as.integer(args[3]) else 300

stopifnot(dir.exists(MODEL_DIR))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

files <- list.files(MODEL_DIR, pattern = "\\.rds$", ignore.case = TRUE, full.names = TRUE)
cat(sprintf("%d fichiers .rds trouves dans %s\n", length(files), MODEL_DIR))

# --- helpers ---------------------------------------------------------------

find_gamlss <- function(x) {
  if (inherits(x, "gamlss")) return(x)
  if (is.list(x)) for (el in x) { g <- find_gamlss(el); if (!is.null(g)) return(g) }
  NULL
}

get_train_data <- function(m, fallback = NULL) {
  d <- tryCatch(eval(m$call$data, envir = environment(formula(m))), error = function(e) NULL)
  if (is.null(d)) d <- tryCatch(eval(m$call$data, envir = globalenv()), error = function(e) NULL)
  if (is.null(d)) d <- fallback
  d
}

# Nom des colonnes de covariables : ADAPTER ICI si different dans tes modeles.
COL_AGE  <- "age"
COL_SEX  <- "sex"
COL_SITE <- "site"
COL_QC   <- "quality"      # ratio de segments rejetes (Autoreject) ; NULL si absent

manifest <- list()
out_rows <- list()

for (f in files) {
  marker <- tools::file_path_sans_ext(basename(f))
  obj <- tryCatch(readRDS(f), error = function(e) NULL)
  m   <- if (is.null(obj)) NULL else find_gamlss(obj)
  if (is.null(m)) { cat(sprintf("[SKIP] %s : pas d'objet gamlss\n", marker)); next }

  dat <- get_train_data(m, fallback = if (is.list(obj)) obj$data else NULL)
  if (is.null(dat)) { cat(sprintf("[SKIP] %s : donnees d'entrainement introuvables\n", marker)); next }

  covs <- intersect(c(COL_AGE, COL_SEX, COL_SITE, COL_QC), names(dat))
  if (!(COL_AGE %in% covs)) { cat(sprintf("[SKIP] %s : colonne age absente\n", marker)); next }

  ages <- dat[[COL_AGE]]
  ages <- ages[is.finite(ages)]
  # grille log-espacee : dense chez les jeunes, ou l'age varie vite
  a_min <- min(ages); a_max <- max(ages)
  shift <- if (a_min <= 0) 1e-3 - a_min else 0
  grid_age <- exp(seq(log(a_min + shift), log(a_max + shift), length.out = N_AGE)) - shift

  sex_levels <- if (COL_SEX %in% covs) sort(unique(as.character(dat[[COL_SEX]]))) else NA_character_
  site_ref   <- if (COL_SITE %in% covs) levels(factor(dat[[COL_SITE]]))[1] else NA_character_
  qc_ref     <- if (COL_QC %in% covs) median(dat[[COL_QC]], na.rm = TRUE) else NA_real_

  nd <- expand.grid(age = grid_age, sex = sex_levels, stringsAsFactors = FALSE)
  names(nd)[1:2] <- c(COL_AGE, COL_SEX)
  if (all(is.na(sex_levels))) nd[[COL_SEX]] <- NULL
  if (COL_SITE %in% covs) nd[[COL_SITE]] <- factor(site_ref, levels = levels(factor(dat[[COL_SITE]])))
  if (COL_QC   %in% covs) nd[[COL_QC]]   <- qc_ref
  if (COL_SEX  %in% covs) nd[[COL_SEX]]  <- factor(nd[[COL_SEX]], levels = levels(factor(dat[[COL_SEX]])))

  p <- tryCatch(
    predictAll(m, newdata = nd, data = dat, type = "response", output = "list"),
    error = function(e) { cat(sprintf("[FAIL] %s : %s\n", marker, conditionMessage(e))); NULL }
  )
  if (is.null(p)) next

  row <- data.frame(
    marker = marker,
    age    = nd[[COL_AGE]],
    sex    = if (COL_SEX %in% covs) as.character(nd[[COL_SEX]]) else NA_character_,
    mu     = as.numeric(p$mu),
    sigma  = if (!is.null(p$sigma)) as.numeric(p$sigma) else NA_real_,
    nu     = if (!is.null(p$nu))    as.numeric(p$nu)    else NA_real_,
    tau    = if (!is.null(p$tau))   as.numeric(p$tau)   else NA_real_,
    stringsAsFactors = FALSE
  )
  out_rows[[marker]] <- row

  # offsets de site : indispensables pour corriger un nouveau sujet enregistre ailleurs
  site_off <- NA
  if (COL_SITE %in% covs) {
    cf <- coef(m, what = "mu")
    site_off <- cf[grepl(paste0("^", COL_SITE), names(cf))]
  }

  manifest[[marker]] <- data.frame(
    marker      = marker,
    family      = m$family[1],
    n_train     = length(m$y),
    age_min     = a_min,
    age_max     = a_max,
    site_ref    = site_ref,
    qc_ref      = qc_ref,
    sex_levels  = paste(sex_levels, collapse = "|"),
    site_offsets = paste(sprintf("%s=%.6g", names(site_off), as.numeric(site_off)), collapse = "|"),
    stringsAsFactors = FALSE
  )
  cat(sprintf("[OK] %s (%s)\n", marker, m$family[1]))
}

params <- do.call(rbind, out_rows)
manif  <- do.call(rbind, manifest)

write.csv(manif, file.path(OUT_DIR, "model_manifest.csv"), row.names = FALSE)
if (requireNamespace("arrow", quietly = TRUE)) {
  arrow::write_parquet(params, file.path(OUT_DIR, "normative_params.parquet"))
} else {
  gz <- gzfile(file.path(OUT_DIR, "normative_params.csv.gz"), "w")
  write.csv(params, gz, row.names = FALSE); close(gz)
}
cat(sprintf("\nTermine : %d modeles exportes -> %s\n", nrow(manif), OUT_DIR))
