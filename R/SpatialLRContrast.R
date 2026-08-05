# ==============================================================================
# SPATIAL LIGAND-RECEPTOR SIGNALLING: BETWEEN-CONDITION CONTRAST
# ------------------------------------------------------------------------------
# Region-of-interest (ROI) aware, patient- and batch-adjusted differential
# analysis of proximity-weighted ligand-receptor signalling. Generalised from a
# glomerulus-level kidney TMA pipeline to ANY user-defined region label.
# ==============================================================================

# ------------------------------------------------------------------------------
# Internal helpers (not exported)
# ------------------------------------------------------------------------------

#' Null-default operator
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Format a number of seconds as a compact human-readable string
#' @noRd
.lr_fmt_time <- function(s) {
  if (!is.finite(s)) return("?")
  if (s < 90)      return(sprintf("%.0fs", s))
  if (s < 5400)    return(sprintf("%.1fm", s / 60))
  sprintf("%.1fh", s / 3600)
}

#' Build a cells x types sparse indicator matrix, keeping ALL requested levels
#'
#' Base \code{model.matrix}/\code{sparse.model.matrix} silently drop factor levels
#' that are absent from the data, which would make the type x type score/contact
#' matrices ragged across regions. This constructs the indicator directly so every
#' region yields the same, full \code{types} columns (empty ones are all-zero).
#' @noRd
.lr_indicator <- function(x, levels) {
  f <- factor(as.character(x), levels = levels)
  j <- as.integer(f)
  keep <- !is.na(j)
  Matrix::sparseMatrix(
    i = which(keep), j = j[keep], x = 1,
    dims = c(length(f), length(levels)),
    dimnames = list(NULL, levels)
  )
}

#' Expand a complex/subunit database into a resolver function
#'
#' Accepts either a data.frame whose row names are complex names and whose cells
#' hold subunit gene symbols (CellChat's \code{$complex} layout; blanks/NA ignored)
#' or a named list of character vectors. Returns a function mapping a character
#' vector of ligand/receptor tokens to their constituent gene symbols.
#' @noRd
.lr_make_resolver <- function(complex_db) {
  if (is.null(complex_db)) return(function(u) unique(u))
  if (is.data.frame(complex_db) || is.matrix(complex_db)) {
    cx <- as.data.frame(complex_db, stringsAsFactors = FALSE)
    lut <- stats::setNames(
      lapply(seq_len(nrow(cx)), function(i) {
        v <- unlist(cx[i, ], use.names = FALSE)
        v <- v[!is.na(v) & nzchar(trimws(as.character(v)))]
        unique(as.character(v))
      }),
      rownames(cx)
    )
  } else if (is.list(complex_db)) {
    lut <- lapply(complex_db, function(v) unique(as.character(v[!is.na(v) & nzchar(v)])))
  } else {
    stop("`complex_db` must be a data.frame, matrix, or named list.")
  }
  function(u) unique(unlist(lapply(u, function(x) if (!is.null(lut[[x]])) lut[[x]] else x),
                            use.names = FALSE))
}

#' Fit one mixed model for a single (interaction, sender, receiver) x contrast
#' @noRd
.lr_fit_one <- function(dat, disease, ref, blocks_shared,
                        has_block, has_patient, use_offset, family_pref) {

  d <- dat[dat$group %in% c(disease, ref), , drop = FALSE]
  if (has_block && !is.null(blocks_shared)) {
    d <- d[as.character(d$block) %in% blocks_shared, , drop = FALSE]
  }
  d <- d[is.finite(d$contacts) & d$contacts > 0, , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d$group <- factor(as.character(d$group), levels = c(ref, disease))
  d <- droplevels(d)
  if (nlevels(d$group) < 2L) return(NULL)

  # per-arm support (used only for reporting counts)
  arms <- split(d, d$group)
  ng   <- vapply(arms, nrow, integer(1))

  multi_block   <- has_block   && length(unique(as.character(d$block)))   > 1L
  multi_patient <- has_patient && length(unique(as.character(d$patient))) > 1L
  d$logc <- log(d$contacts)

  fixed <- "group"
  if (multi_block) fixed <- paste(fixed, "+ block")
  reff  <- if (multi_patient) " + (1 | patient)" else ""
  offs  <- if (use_offset) " + offset(logc)" else ""

  fam <- NA_character_; m <- NULL
  if (identical(family_pref, "tweedie")) {
    form <- stats::as.formula(paste0("score ~ ", fixed, reff, offs))
    # suppressWarnings (not a warning handler): a convergence warning should not
    # abort the fit -- doing so also trips glmmTMB's internal system.time and
    # prints "Timing stopped at:". We judge convergence from pdHess below.
    m <- tryCatch(suppressWarnings(glmmTMB::glmmTMB(form, data = d,
                                   family = glmmTMB::tweedie(link = "log"))),
                  error = function(e) NULL)
    fam <- "tweedie"
    ok <- !is.null(m) && isTRUE(tryCatch(m$sdr$pdHess, error = function(e) FALSE))
    if (!ok) m <- NULL
  }
  if (is.null(m)) {   # Gaussian(log1p) fallback (or primary if family_pref == "gaussian")
    d$y <- if (use_offset) log1p(d$score / d$contacts) else log1p(d$score)
    form2 <- stats::as.formula(paste0("y ~ ", fixed, reff))
    m <- tryCatch(suppressWarnings(glmmTMB::glmmTMB(form2, data = d)),
                  error = function(e) NULL)
    fam <- "gaussian_log1p"
    if (is.null(m)) return(NULL)
  }

  co <- broom.mixed::tidy(m, effects = "fixed")
  co <- co[startsWith(co$term, paste0("group", disease)), , drop = FALSE]
  if (!nrow(co)) return(NULL)

  data.frame(
    disease = disease, ref = ref, family = fam,
    estimate = co$estimate[1], std.error = co$std.error[1], p = co$p.value[1],
    n_glom = nrow(d),
    n_pat  = if (has_patient) length(unique(d$patient)) else NA_integer_,
    n_glom_disease = as.integer(ng[disease]),
    n_glom_ref     = as.integer(ng[ref]),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# Main function
# ------------------------------------------------------------------------------

#' Between-condition contrast of spatial ligand-receptor signalling in regions of interest
#'
#' Tests whether proximity-weighted ligand-receptor (L-R) signalling between two
#' cell types \emph{differs between experimental conditions}, using the region of
#' interest (ROI) as the unit of observation and accounting for repeated sampling
#' of patients and for batch (e.g. tissue-microarray / slide) structure. It is a
#' between-group companion to within-group L-R enrichment tests: rather than asking
#' "is sender\eqn{\rightarrow}receiver signalling enriched against a within-region
#' null", it asks "does per-contact sender\eqn{\rightarrow}receiver signalling
#' \emph{change} between conditions".
#'
#' Although developed on kidney glomeruli in tissue microarrays, the region is
#' fully generic: any categorical cell-level label that partitions cells into
#' spatial units (glomeruli, tumour nests, follicles, crypts, niches, manually
#' drawn AOIs, ...) can be supplied via \code{region_col}. All variables are set
#' through arguments, so the same routine covers other designs.
#'
#' @details
#' \strong{Pipeline.}
#' \enumerate{
#'   \item \emph{Proximity graph.} Within each region, cells closer than
#'     \code{radius} (in the units of \code{coord_cols} after multiplying by
#'     \code{coord_scale}) are joined by an unweighted, undirected edge. Regions
#'     with fewer than \code{min_cells_region} cells, or with no proximate pair,
#'     are dropped.
#'   \item \emph{Testable L-R set (cell-type-aware prevalence filter).} A gene is
#'     retained if, in \emph{at least one} sufficiently large cell type
#'     (\eqn{\ge}\code{min_type_cells} cells pooled over kept regions), it is
#'     detected in \eqn{\ge}\code{prev_min} of that type's cells AND in
#'     \eqn{\ge}\code{min_expr_cells} cells. This rescues cell-type-restricted
#'     ligands/receptors that a global filter dilutes away, while still producing a
#'     \emph{single} tested set applied identically to every condition (so the
#'     between-condition contrasts remain valid). Set \code{prevalence = "global"}
#'     for a pooled filter or \code{"none"} to skip filtering.
#'   \item \emph{Per-region scores and contacts.} For each L-R pair and each
#'     ordered (sender, receiver) cell-type combination the function computes an
#'     interaction \code{score} and the number of proximal sender-receiver
#'     \code{contacts} (used as the model offset). With \code{scoring = "sum"}
#'     (default, CellChat-like additive score)
#'     \deqn{s_{a\to b}=\sum_{i\in a}L_i\,\#\{b\text{-neighbours of }i\}+\sum_{j\in b}R_j\,\#\{a\text{-neighbours of }j\},}
#'     with \eqn{L,R} the mean expression across ligand/receptor subunits. With
#'     \code{scoring = "product"}, \eqn{s_{a\to b}=\sum_{i\in a,\,j\in b,\,i\sim j}L_i R_j}.
#'   \item \emph{Model.} For every (interaction, sender, receiver) \eqn{\times}
#'     contrast a generalized linear mixed model is fit:
#'     \code{score ~ group + block + (1 | patient) + offset(log(contacts))}. The
#'     \code{group} coefficient is the condition effect; with the offset it is a
#'     \emph{per-contact} intensity contrast, so a condition is not called "more
#'     signalling" merely for containing more of a cell type (drop the offset with
#'     \code{use_offset = FALSE} to interpret it as \emph{total} signalling). The
#'     block term and patient random intercept are included only when they vary in
#'     the contrast's data. The primary family is Tweedie(log) (non-negative,
#'     continuous, exact zeros); on non-convergence the fit falls back to a
#'     Gaussian LMM on \code{log1p} of the per-contact score. Multiplicity is
#'     controlled within each contrast by \code{padj_method}.
#' }
#'
#' \strong{Batch (block) handling.} When batch is confounded with condition
#' (e.g. each TMA carries only some diagnoses), set \code{match_blocks = TRUE}
#' (default): each contrast is restricted to the blocks shared by its two groups,
#' turning batch into a balanced within-contrast fixed block rather than a
#' confound. Contrasts whose two groups share \emph{no} block are skipped with a
#' message. With a balanced design (batch not confounded with condition) set
#' \code{match_blocks = FALSE} to use all data.
#'
#' @param object A \code{Seurat} object. Expression is taken from
#'   \code{GetAssayData(object, assay, layer)}; per-cell metadata from
#'   \code{object[[]]}. May be \code{NULL} if \code{expr}, \code{meta} (and
#'   \code{coords}) are supplied directly (non-Seurat use).
#' @param expr Optional genes \eqn{\times} cells matrix (dense or
#'   \pkg{Matrix}-sparse) overriding extraction from \code{object}. Column order
#'   must match \code{meta}/\code{coords}.
#' @param meta Optional per-cell metadata \code{data.frame} (one row per cell,
#'   same order as \code{expr}) overriding \code{object[[]]}.
#' @param coords Optional cells \eqn{\times} 2 numeric matrix/data.frame of
#'   coordinates overriding \code{coord_cols}. If \code{NULL}, coordinates are read
#'   from \code{coord_cols} in the metadata, or, failing that, from
#'   \code{SeuratObject::GetTissueCoordinates(object)}.
#' @param region_col Metadata column giving each cell's region/ROI label
#'   (the generalisation of \code{glom_id}). Cells with \code{NA} here are ignored.
#' @param celltype_col Metadata column giving each cell's type/cluster.
#' @param group_col Metadata column giving the experimental condition to contrast.
#' @param patient_col Metadata column identifying the biological replicate
#'   (patient/donor) used as a random intercept. \code{NULL} drops the random
#'   effect (fixed-effects only); a warning notes the pseudo-replication risk.
#' @param block_col Metadata column giving the batch/block (e.g. TMA or slide).
#'   \code{NULL} disables block matching and the block fixed effect.
#' @param assay,layer Seurat assay and layer/slot to pull expression from.
#'   Defaults \code{"RNA"} and \code{"counts"}.
#' @param coord_cols Length-2 character vector naming the x/y metadata columns.
#'   Default \code{c("x", "y")}. Ignored if \code{coords} is supplied.
#' @param coord_scale Multiplier converting coordinate units so that
#'   \code{radius} is expressed in the same units (e.g. set to the microns-per-unit
#'   factor if coordinates are not already in microns). Default 1.
#' @param lr_db Ligand-receptor database. Either a \code{data.frame} with the
#'   columns named in \code{lr_cols}, or a list carrying such a frame in
#'   \code{$interaction} (and, optionally, a complex table in \code{$complex}), as
#'   returned by CellChat (\code{CellChatDB.human}). Required.
#' @param complex_db Optional complex/subunit table resolving multi-subunit
#'   ligands/receptors: a \code{data.frame} whose row names are complex names and
#'   whose cells hold subunit symbols (CellChat's \code{$complex}), or a named list
#'   of subunit vectors. If \code{lr_db} is a CellChat-style list its
#'   \code{$complex} is used automatically unless \code{complex_db} is given.
#' @param lr_cols Named character vector mapping the roles
#'   \code{interaction, ligand, receptor, pathway} to column names in the L-R
#'   frame. Default matches CellChat:
#'   \code{c(interaction = "interaction_name", ligand = "ligand",
#'   receptor = "receptor", pathway = "pathway_name")}. If a \code{pathway} column
#'   is absent it is set to \code{NA}.
#' @param radius Proximity distance defining a contact, in coordinate units after
#'   \code{coord_scale}. Default 25.
#' @param min_cells_region Minimum cells for a region to be analysed. Default 20.
#' @param prevalence One of \code{"celltype"} (default, per-cell-type filter),
#'   \code{"global"} (pooled fraction over all kept region cells), or \code{"none"}.
#' @param prev_min Minimum detection fraction for the prevalence filter. Default
#'   0.05. Re-run at 0.10 as a sensitivity check.
#' @param min_expr_cells Minimum number of expressing cells for the prevalence
#'   filter. Default 10.
#' @param min_type_cells Minimum cells for a cell type to be trusted by the
#'   \code{"celltype"} filter. Default 30.
#' @param scoring \code{"sum"} (default, additive) or \code{"product"}
#'   (pairwise \eqn{L_i R_j}).
#' @param use_offset If \code{TRUE} (default) include \code{offset(log(contacts))}
#'   for a per-contact intensity contrast; if \code{FALSE} omit it for a total
#'   signalling contrast.
#' @param reference_group Condition treated as the baseline when \code{contrasts}
#'   is \code{NULL}. Defaults to \code{"Normal"} if present, else the first level
#'   of \code{group_col}.
#' @param contrasts Optional named list of contrasts; each element is a list with
#'   \code{disease} and optional \code{ref} (defaults to \code{reference_group}),
#'   e.g. \code{list("MN_vs_Normal" = list(disease = "MN"),
#'   "FSGS_vs_MCD" = list(disease = "FSGS", ref = "MCD"))}. If \code{NULL}, every
#'   non-reference group is contrasted against \code{reference_group}.
#' @param match_blocks If \code{TRUE} (default) restrict each contrast to the
#'   blocks shared by its two groups (see Details). Ignored when \code{block_col}
#'   is \code{NULL}.
#' @param min_region Minimum informative regions per arm to attempt a fit.
#'   Default 4.
#' @param min_patient Minimum patients per arm to attempt a fit. Default 2
#'   (forced to 1 when \code{patient_col} is \code{NULL}).
#' @param family \code{"tweedie"} (default; Gaussian-log1p fallback) or
#'   \code{"gaussian"} (Gaussian-log1p only).
#' @param padj_method Multiple-testing adjustment passed to \code{stats::p.adjust}
#'   \emph{within each contrast}. Default \code{"BH"}.
#' @param parallel If \code{TRUE}, parallelise the per-pair fits with
#'   \code{parallel::mclapply} (forking; not on Windows). Default \code{FALSE}.
#' @param n_cores Number of cores when \code{parallel = TRUE}. Default 1.
#' @param outdir Optional directory. If non-\code{NULL}, the tested pair table,
#'   per-region scores, results table and (if \code{make_volcano}) a volcano PDF
#'   are written there. If \code{NULL}, nothing is written and results are only
#'   returned.
#' @param make_volcano If \code{TRUE} (default) and \pkg{ggplot2} is available,
#'   build a per-contrast volcano plot (returned, and saved when \code{outdir} is
#'   set).
#' @param seed RNG seed for reproducibility of any stochastic model internals.
#'   Default 1.
#' @param verbose Print progress messages. Default \code{TRUE}.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{results}}{tidy \code{data.frame}, one row per tested
#'       (contrast, interaction, sender, receiver): \code{estimate},
#'       \code{std.error}, \code{p}, \code{padj}, \code{log2FC}
#'       (\eqn{= \mathrm{estimate}/\log 2} under the log link), family used, and
#'       per-arm support counts. Sorted by contrast, \code{padj}, |\code{log2FC}|.}
#'     \item{\code{scores}}{per-region \code{data.frame} of \code{score} and
#'       \code{contacts} with region metadata.}
#'     \item{\code{tested_lr}}{the tested L-R set with resolved ligand/receptor
#'       genes.}
#'     \item{\code{volcano}}{the \pkg{ggplot} object (or \code{NULL}).}
#'     \item{\code{params}}{the resolved arguments, for the methods section.}
#'   }
#'
#' @section Reproducibility / reviewer notes:
#' The observation unit is the region and patient is a random intercept, so the
#' test is not driven by one patient contributing many regions. The
#' \code{log(contacts)} offset makes the contrast per-contact. Record which
#' \code{family} was used per row. Recommended sensitivity checks: vary
#' \code{radius} (e.g. 15 and 35), compare \code{scoring} and \code{use_offset},
#' re-run at \code{prev_min = 0.10}, and confirm top hits with an orthogonal
#' readout.
#'
#' @seealso \code{\link{AssignCellsToAOIs}}, \code{\link{ReadQuPathROIs}} for
#'   producing the region label consumed by \code{region_col}.
#'
#' @references
#' Jin S. et al. Inference and analysis of cell-cell communication using CellChat.
#' \emph{Nature Communications} (2021). Brooks M.E. et al. glmmTMB balances speed
#' and flexibility among packages for zero-inflated generalized linear mixed
#' models. \emph{The R Journal} (2017).
#'
#' @examples
#' \dontrun{
#' library(gbspatial)
#' data("my_spatial_data")
#'
#' # CellChat provides both the interaction table and the complex table:
#' # CellChatDB.human$interaction / $complex
#' res <- SpatialLRContrast(
#'   object       = my_spatial_data,
#'   region_col   = "aoi_id",           # any ROI label; was glom_id
#'   celltype_col = "cell_type_high_res",
#'   group_col    = "group",
#'   patient_col  = "id",
#'   block_col    = "slide",            # TMA / batch
#'   assay = "RNA", layer = "counts",
#'   coord_cols = c("x", "y"),
#'   lr_db = CellChatDB.human,          # list with $interaction + $complex
#'   radius = 25,
#'   contrasts = list(
#'     "MN_vs_Normal"   = list(disease = "MN"),
#'     "FSGS_vs_Normal" = list(disease = "FSGS-tip"),
#'     "FSGS_vs_MCD"    = list(disease = "FSGS-tip", ref = "MCD")
#'   ),
#'   outdir = tempdir()
#' )
#' head(res$results)
#' res$volcano
#' }
#'
#' @importFrom stats as.formula dist p.adjust setNames
#' @importFrom utils write.csv
#' @export
SpatialLRContrast <- function(object = NULL,
                              expr = NULL, meta = NULL, coords = NULL,
                              region_col, celltype_col, group_col,
                              patient_col = NULL, block_col = NULL,
                              assay = "RNA", layer = "counts",
                              coord_cols = c("x", "y"), coord_scale = 1,
                              lr_db,
                              complex_db = NULL,
                              lr_cols = c(interaction = "interaction_name",
                                          ligand = "ligand",
                                          receptor = "receptor",
                                          pathway = "pathway_name"),
                              radius = 25, min_cells_region = 20,
                              prevalence = c("celltype", "global", "none"),
                              prev_min = 0.05, min_expr_cells = 10,
                              min_type_cells = 30,
                              scoring = c("sum", "product"),
                              use_offset = TRUE,
                              reference_group = NULL,
                              contrasts = NULL, match_blocks = TRUE,
                              min_region = 4, min_patient = 2,
                              family = c("tweedie", "gaussian"),
                              padj_method = "BH",
                              parallel = FALSE, n_cores = 1,
                              outdir = NULL, make_volcano = TRUE,
                              seed = 1, verbose = TRUE) {

  prevalence <- match.arg(prevalence)
  scoring    <- match.arg(scoring)
  family     <- match.arg(family)
  for (p in c("Matrix", "glmmTMB", "broom.mixed"))
    if (!requireNamespace(p, quietly = TRUE))
      stop("Package '", p, "' is required for SpatialLRContrast().")
  set.seed(seed)
  say <- function(...) if (isTRUE(verbose)) message(...)

  # ---- resolve expression, metadata, coordinates --------------------------
  if (is.null(expr) || is.null(meta)) {
    if (is.null(object)) stop("Supply either `object`, or `expr` and `meta` (and `coords`).")
    if (!requireNamespace("SeuratObject", quietly = TRUE))
      stop("Package 'SeuratObject' is required to read from a Seurat object.")
    if (is.null(expr)) expr <- SeuratObject::GetAssayData(object, assay = assay, layer = layer)
    if (is.null(meta)) meta <- as.data.frame(object[[]])
  }
  meta <- as.data.frame(meta)
  if (ncol(expr) != nrow(meta))
    stop("`expr` has ", ncol(expr), " columns but `meta` has ", nrow(meta), " rows.")

  need <- c(region_col, celltype_col, group_col)
  if (!is.null(patient_col)) need <- c(need, patient_col)
  if (!is.null(block_col))   need <- c(need, block_col)
  miss <- setdiff(need, colnames(meta))
  if (length(miss)) stop("Missing metadata column(s): ", paste(miss, collapse = ", "))

  if (is.null(coords)) {
    if (all(coord_cols %in% colnames(meta))) {
      coords <- as.matrix(meta[, coord_cols, drop = FALSE])
    } else if (!is.null(object) && requireNamespace("SeuratObject", quietly = TRUE)) {
      tc <- SeuratObject::GetTissueCoordinates(object)
      xy <- intersect(c("x", "y"), colnames(tc))
      if (length(xy) < 2) stop("Could not find x/y in GetTissueCoordinates(); set `coord_cols` or pass `coords`.")
      coords <- as.matrix(tc[, xy, drop = FALSE])
      if (nrow(coords) != nrow(meta))
        stop("Tissue coordinates (", nrow(coords), ") do not align with cells (", nrow(meta), ").")
    } else {
      stop("Coordinates not found. Provide `coords`, or `coord_cols` present in metadata.")
    }
  } else {
    coords <- as.matrix(coords)
  }
  storage.mode(coords) <- "double"
  coords <- coords * coord_scale
  if (nrow(coords) != nrow(meta)) stop("`coords` rows must equal number of cells.")

  cell_type <- as.character(meta[[celltype_col]])
  region    <- as.character(meta[[region_col]])
  types_all <- sort(unique(cell_type[!is.na(cell_type)]))

  if (is.null(patient_col)) {
    warning("`patient_col` is NULL: fitting without a patient random effect. ",
            "Regions from the same donor are treated as independent, which risks pseudo-replication.")
    has_patient <- FALSE
    meta$.patient <- factor("all")
    min_patient <- 1L
  } else {
    has_patient <- TRUE
    meta$.patient <- as.character(meta[[patient_col]])
  }
  has_block <- !is.null(block_col)
  meta$.block <- if (has_block) as.character(meta[[block_col]]) else "all"

  # ---- regions and proximity graphs ---------------------------------------
  cell_idx  <- which(!is.na(region))
  region_list <- split(cell_idx, region[cell_idx])

  glom_W <- lapply(region_list, function(idx) {
    if (length(idx) < min_cells_region) return(NULL)
    d <- as.matrix(stats::dist(coords[idx, , drop = FALSE]))
    W <- (d > 0 & d <= radius) * 1
    if (sum(W) == 0) return(NULL)
    Matrix::Matrix(W, sparse = TRUE)
  })
  keep <- !vapply(glom_W, is.null, logical(1))
  glom_W      <- glom_W[keep]
  valid_gloms <- names(glom_W)
  say(sprintf("Usable regions (>=%d cells, >=1 proximate pair): %d of %d",
              min_cells_region, length(valid_gloms), length(region_list)))
  if (!length(valid_gloms)) stop("No region passed the size / proximity filter.")

  # ---- region-level metadata ----------------------------------------------
  first_by_region <- function(v) vapply(region_list[valid_gloms],
                                        function(idx) as.character(v[idx[1]]), character(1))
  gmeta <- data.frame(
    glom_id = valid_gloms,
    group   = factor(first_by_region(meta[[group_col]])),
    patient = first_by_region(meta$.patient),
    block   = first_by_region(meta$.block),
    stringsAsFactors = FALSE
  )

  # ---- L-R database ---------------------------------------------------------
  if (is.list(lr_db) && !is.data.frame(lr_db) && !is.null(lr_db$interaction)) {
    if (is.null(complex_db) && !is.null(lr_db$complex)) complex_db <- lr_db$complex
    lr_tab <- lr_db$interaction
  } else {
    lr_tab <- lr_db
  }
  lr_tab <- as.data.frame(lr_tab, stringsAsFactors = FALSE)
  for (r in c("interaction", "ligand", "receptor"))
    if (!lr_cols[[r]] %in% colnames(lr_tab))
      stop("L-R table lacks the '", r, "' column '", lr_cols[[r]], "' (set via `lr_cols`).")
  has_path <- !is.null(lr_cols[["pathway"]]) && lr_cols[["pathway"]] %in% colnames(lr_tab)
  resolve  <- .lr_make_resolver(complex_db)

  # ---- prevalence filter ----------------------------------------------------
  glom_cells <- unlist(region_list[valid_gloms], use.names = FALSE)
  gene_names <- rownames(expr)
  if (identical(prevalence, "none")) {
    keepgene <- function(gs) intersect(gs, gene_names)
  } else if (identical(prevalence, "global")) {
    pos <- (expr[, glom_cells, drop = FALSE] > 0)
    n_expr <- Matrix::rowSums(pos)
    frac   <- n_expr / length(glom_cells)
    names(n_expr) <- names(frac) <- gene_names
    keepgene <- function(gs) {
      gs <- intersect(gs, gene_names); if (!length(gs)) return(character(0))
      gs[frac[gs] >= prev_min & n_expr[gs] >= min_expr_cells]
    }
  } else { # celltype-aware
    ctv <- .lr_indicator(cell_type[glom_cells], types_all)     # cells x types
    n_ct      <- Matrix::colSums(ctv)
    pos       <- (expr[, glom_cells, drop = FALSE] > 0) * 1
    n_expr_ct <- as.matrix(pos %*% ctv)                        # genes x type
    detfrac   <- sweep(n_expr_ct, 2, pmax(n_ct, 1), "/")
    rownames(detfrac) <- rownames(n_expr_ct) <- gene_names
    ok_types  <- names(n_ct)[n_ct >= min_type_cells]
    if (!length(ok_types))
      stop("No cell type reaches `min_type_cells` (", min_type_cells, ").")
    keepgene <- function(gs) {
      gs <- intersect(gs, gene_names); if (!length(gs)) return(character(0))
      f <- detfrac[gs, ok_types, drop = FALSE]
      n <- n_expr_ct[gs, ok_types, drop = FALSE]
      gs[rowSums(f >= prev_min & n >= min_expr_cells) > 0]
    }
  }

  spl <- function(x) unique(unlist(strsplit(as.character(x), "[,_&|+/ ]+")))
  lig_raw <- lapply(lr_tab[[lr_cols[["ligand"]]]],   function(z) resolve(spl(z)))
  rec_raw <- lapply(lr_tab[[lr_cols[["receptor"]]]], function(z) resolve(spl(z)))
  lig <- lapply(lig_raw, keepgene)
  rec <- lapply(rec_raw, keepgene)
  ok  <- lengths(lig) > 0 & lengths(rec) > 0
  lr_use <- data.frame(
    interaction = as.character(lr_tab[[lr_cols[["interaction"]]]])[ok],
    pathway     = if (has_path) as.character(lr_tab[[lr_cols[["pathway"]]]])[ok] else NA_character_,
    stringsAsFactors = FALSE
  )
  lr_use$lig <- lig[ok]; lr_use$rec <- rec[ok]
  lr_use <- lr_use[!duplicated(lr_use$interaction), , drop = FALSE]
  if (!nrow(lr_use)) stop("No L-R pair survived the prevalence filter; loosen `prev_min`/`min_expr_cells`.")
  say(nrow(lr_use), " L-R pairs testable (prevalence = ", prevalence,
      "; prev_min = ", prev_min, ", min_expr_cells = ", min_expr_cells, ")")

  tested_lr <- data.frame(
    interaction = lr_use$interaction, pathway = lr_use$pathway,
    ligand_genes   = vapply(lr_use$lig, paste, character(1), collapse = ";"),
    receptor_genes = vapply(lr_use$rec, paste, character(1), collapse = ";"),
    stringsAsFactors = FALSE
  )

  # ---- per-region scores + contacts ----------------------------------------
  perglom_one <- function(gm) {
    idx <- region_list[[gm]]
    if (length(idx) < 3L) return(NULL)
    d <- as.matrix(stats::dist(coords[idx, , drop = FALSE]))
    W <- (d > 0 & d <= radius) * 1
    if (sum(W) == 0) return(NULL)
    P  <- .lr_indicator(cell_type[idx], types_all)      # cells x type
    WP <- as.matrix(W %*% P)                            # cells x type: neighbour counts
    contacts <- as.matrix(Matrix::crossprod(P, WP))     # type x type: proximal pairs
    ij <- which(contacts > 0, arr.ind = TRUE)
    if (!nrow(ij)) return(NULL)

    e <- as.matrix(expr[unique(unlist(c(lr_use$lig, lr_use$rec))), idx, drop = FALSE])
    out <- vector("list", nrow(lr_use))
    for (i in seq_len(nrow(lr_use))) {
      L <- colMeans(e[lr_use$lig[[i]], , drop = FALSE])
      R <- colMeans(e[lr_use$rec[[i]], , drop = FALSE])
      if (identical(scoring, "sum")) {
        score <- as.matrix(Matrix::crossprod(P, L * WP) + Matrix::crossprod(WP, R * P))
      } else { # product
        A <- Matrix::Diagonal(x = L) %*% P
        B <- Matrix::Diagonal(x = R) %*% P
        score <- as.matrix(Matrix::crossprod(A, W %*% B))
      }
      out[[i]] <- data.frame(
        glom_id = gm, interaction = lr_use$interaction[i], pathway = lr_use$pathway[i],
        sender = types_all[ij[, 1]], receiver = types_all[ij[, 2]],
        score = score[ij], contacts = contacts[ij], stringsAsFactors = FALSE)
    }
    do.call(rbind, out)
  }

  lapply_fun <- if (isTRUE(parallel) && requireNamespace("parallel", quietly = TRUE))
    function(X, FUN) parallel::mclapply(X, FUN, mc.cores = n_cores) else lapply
  perglom <- do.call(rbind, lapply_fun(valid_gloms, perglom_one))
  perglom <- merge(perglom, gmeta, by = "glom_id", all.x = TRUE)
  say("per-region rows: ", nrow(perglom))

  # ---- contrasts ------------------------------------------------------------
  glevels <- levels(gmeta$group)
  if (is.null(reference_group))
    reference_group <- if ("Normal" %in% glevels) "Normal" else glevels[1]
  if (is.null(contrasts)) {
    others <- setdiff(glevels, reference_group)
    contrasts <- stats::setNames(
      lapply(others, function(g) list(disease = g, ref = reference_group)),
      paste0(others, "_vs_", reference_group))
  }
  blocks_of <- function(grp) unique(as.character(gmeta$block[gmeta$group == grp]))

  # Index perglom ONCE by (interaction, sender, receiver) so each tested
  # combination is fetched directly instead of re-scanning all ~millions of
  # per-region rows for every combination (the previous O(rows x combinations)
  # cost was the dominant runtime).
  key        <- paste(perglom$interaction, perglom$sender, perglom$receiver, sep = "\u001f")
  row_groups <- split(seq_len(nrow(perglom)), key)
  first_idx  <- vapply(row_groups, function(z) z[1L], integer(1))
  grp_meta   <- perglom[first_idx, c("interaction", "pathway", "sender", "receiver"), drop = FALSE]
  rownames(grp_meta) <- NULL
  n_groups   <- length(row_groups)
  say(sprintf("Testing up to %d cell-type pair(s) per contrast across %d contrast(s)%s.",
              n_groups, length(contrasts),
              if (isTRUE(parallel)) sprintf(" (parallel: %d cores)", n_cores) else ""))

  # fit one combination k for the current contrast (returns a 1-row df or NULL)
  fit_group <- function(k, disease, ref, blocks_shared, cname) {
    dat <- perglom[row_groups[[k]], , drop = FALSE]
    dd  <- dat[dat$group %in% c(disease, ref) & dat$contacts > 0, , drop = FALSE]
    if (has_block && isTRUE(match_blocks))
      dd <- dd[as.character(dd$block) %in% blocks_shared, , drop = FALSE]
    if (!nrow(dd)) return(NULL)
    gf <- factor(as.character(dd$group), levels = c(ref, disease))
    ng <- tabulate(gf, nbins = 2L); names(ng) <- c(ref, disease)
    if (any(ng < min_region)) return(NULL)
    if (has_patient) {
      np <- c(length(unique(dd$patient[gf == ref])), length(unique(dd$patient[gf == disease])))
      if (any(np < min_patient)) return(NULL)
    }
    r <- .lr_fit_one(dat, disease = disease, ref = ref, blocks_shared = blocks_shared,
                     has_block = has_block, has_patient = has_patient,
                     use_offset = use_offset, family_pref = family)
    if (is.null(r)) return(NULL)
    cbind(data.frame(contrast = cname, stringsAsFactors = FALSE),
          grp_meta[k, , drop = FALSE], r)
  }

  results <- list()
  for (cname in names(contrasts)) {
    spec <- contrasts[[cname]]
    disease <- spec$disease; ref <- spec$ref %||% reference_group
    blocks_shared <- if (has_block && isTRUE(match_blocks))
      intersect(blocks_of(disease), blocks_of(ref)) else NULL
    if (has_block && isTRUE(match_blocks) && length(blocks_shared) == 0L) {
      say("== ", cname, " == skipped: groups share no block (batch confounded).")
      next
    }
    say("== ", cname, " ==")

    if (isTRUE(parallel) && requireNamespace("parallel", quietly = TRUE)) {
      res_c <- parallel::mclapply(seq_len(n_groups), fit_group,
                                  disease = disease, ref = ref,
                                  blocks_shared = blocks_shared, cname = cname,
                                  mc.cores = n_cores)
    } else {
      res_c <- vector("list", n_groups)
      t0    <- proc.time()[["elapsed"]]
      step  <- max(1L, n_groups %/% 200L)          # ~200 progress updates
      for (k in seq_len(n_groups)) {
        res_c[[k]] <- fit_group(k, disease, ref, blocks_shared, cname)
        if (isTRUE(verbose) && (k %% step == 0L || k == n_groups)) {
          el   <- proc.time()[["elapsed"]] - t0
          frac <- k / n_groups
          eta  <- if (frac > 0) el * (1 - frac) / frac else NA_real_
          cat(sprintf("\r  %s  %5.1f%%  (%d/%d)  elapsed %s  ETA %s      ",
                      cname, 100 * frac, k, n_groups,
                      .lr_fmt_time(el), .lr_fmt_time(eta)))
          utils::flush.console()
        }
      }
      if (isTRUE(verbose)) cat("\n")
    }
    res_c <- do.call(rbind, res_c)
    if (!is.null(res_c) && nrow(res_c)) res_c$padj <- stats::p.adjust(res_c$p, padj_method)
    results[[cname]] <- res_c
  }
  between <- do.call(rbind, results)
  if (!is.null(between) && nrow(between)) {
    between$log2FC <- between$estimate / log(2)
    between <- between[order(between$contrast, between$padj, -abs(between$log2FC)), ]
    rownames(between) <- NULL
  } else {
    between <- data.frame()
  }

  # ---- volcano --------------------------------------------------------------
  volcano <- NULL
  if (isTRUE(make_volcano) && nrow(between) &&
      requireNamespace("ggplot2", quietly = TRUE)) {
    volcano <- ggplot2::ggplot(between,
        ggplot2::aes(x = .data$log2FC, y = -log10(.data$padj))) +
      ggplot2::geom_point(ggplot2::aes(color = .data$padj < 0.05),
                          size = 0.8, alpha = 0.7) +
      ggplot2::geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
      ggplot2::scale_color_manual(values = c("grey70", "firebrick"), guide = "none") +
      ggplot2::facet_wrap(~ contrast, scales = "free") +
      ggplot2::labs(x = "log2 fold change (disease vs reference, per contact)",
                    y = expression(-log[10]~padj)) +
      ggplot2::theme_bw(9)
  }

  # ---- optional write-out ---------------------------------------------------
  if (!is.null(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(tested_lr, file.path(outdir, "tested_lr_pairs.csv"), row.names = FALSE)
    saveRDS(perglom, file.path(outdir, "perregion_scores.rds"))
    utils::write.csv(between, file.path(outdir, "spatial_lr_between_condition.csv"),
                     row.names = FALSE)
    if (!is.null(volcano))
      ggplot2::ggsave(file.path(outdir, "volcano_between_condition.pdf"),
                      volcano, width = 10, height = 8)
    say("Wrote results to ", outdir)
  }

  if (nrow(between)) {
    sig <- between[between$padj < 0.05 & !is.na(between$padj), , drop = FALSE]
    say("Differential (padj<0.05) per-contact L-R interactions per contrast:")
    if (isTRUE(verbose)) print(table(sig$contrast))
  }

  list(
    results   = between,
    scores    = perglom,
    tested_lr = tested_lr,
    volcano   = volcano,
    params    = list(region_col = region_col, celltype_col = celltype_col,
                     group_col = group_col, patient_col = patient_col,
                     block_col = block_col, assay = assay, layer = layer,
                     radius = radius, coord_scale = coord_scale,
                     min_cells_region = min_cells_region, prevalence = prevalence,
                     prev_min = prev_min, min_expr_cells = min_expr_cells,
                     min_type_cells = min_type_cells, scoring = scoring,
                     use_offset = use_offset, reference_group = reference_group,
                     match_blocks = match_blocks, min_region = min_region,
                     min_patient = min_patient, family = family,
                     padj_method = padj_method, n_regions_used = length(valid_gloms),
                     n_lr_tested = nrow(tested_lr), cell_types = types_all)
  )
}
