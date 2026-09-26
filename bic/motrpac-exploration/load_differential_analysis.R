#' @title Load Differential Analysis Results
#'
#' @param selected_omes character; one of \code{\link{ome_available_list}}.
#' @param selected_tissues character; one of
#'   \code{\link{tissue_available_list}}.
#' @param single_matrix logical; if \code{TRUE}, returns a single
#'   \code{data.frame} containing all results. Otherwise, returns a list of
#'   \code{data.frame} objects (default).
#' @param epigen logical; a toggle of TRUE/FALSE if epigenetics data is desired.
#'   The epigenomics tables are not shipped in the package; they are downloaded
#'   from the public c2.0 release on the MoTrPAC CloudFront distribution at run
#'   time. No credentials are needed, but downloads are not cached and are slow
#'   due to file sizes.
#' @param repo_local_dir Deprecated and ignored. Epigenomics files are no longer
#'   downloaded to a local cache; supplying it only prints a message.
#' @param combine_with_featgene logical; whether to include columns from
#'   \code{HUMAN_FEATURE_TO_GENE} in the output.
#' @param verbose logical; whether or not to display messages for some warnings.
#'
#' @param load_clinical logical; whether to include the clinical chemistry omes
#'   (\code{clinical_ome_list()}: \code{"prot-clinical"} and
#'   \code{"metab-t-clinical"}). \code{FALSE} by default, so \code{"all"}
#'   returns the research omes and nothing changes for callers written before
#'   v2.0 split clinical chemistry out. Set \code{TRUE} to include them; they
#'   are dropped even when named unless it is set.
#'
#' @returns A nested list of \code{data.table} objects. The top level names are
#'   the tissues, while the second level names are the omes. Each table may
#'   possess the following columns:
#'
#'   \describe{
#'     \item{tissue}{factor; the tissue.}
#'     \item{assay}{factor; the assay family, not the platform. Every
#'     metabolomics table reads \code{"metab"} — clinical chemistry included —
#'     so \code{assay} alone does not separate \code{metab-t-clinical} from the
#'     research platforms, and five analytes (Cortisol, Glucose, Glycerol, KET,
#'     NEFA) exist on both. Include \code{platform} in any key that has to tell
#'     them apart.}
#'     \item{platform}{factor; (metabolomics only) the metabolomics platform.}
#'     \item{full_model}{factor; full model containing predictors and any
#'     covariates.}
#'     \item{contrast}{factor; full contrast (up to 33).}
#'     \item{contrast_short}{factor; shortened version of the contrasts.}
#'     \item{contrast_type}{factor; one of "exercise_with_controls",
#'     "exercise_no_controls", "Endur_vs_Resist", "baseline", or
#'     "control_only".}
#'     \item{contrast_category}{factor; one of "EE-CON", "RE-CON", "EE-EE",
#'     "RE-RE", "EE-RE", or "CON-CON".}
#'     \item{feature_id}{factor; feature identifier (may be proteins,
#'     phosphosites, transcripts, metabolites/lipids, peaks, or GpGs).}
#'     \item{logFC}{numeric; difference between the group means in the
#'     contrast.}
#'     \item{CI.L_calculated}{numeric; lower bound of the 95\% confidence
#'     interval on \code{logFC}, computed as
#'     \code{logFC - (logFC / t) * qt(0.975, df)} against that contrast's own
#'     residual degrees of freedom.}
#'     \item{CI.R_calculated}{numeric; upper bound of the same interval. The
#'     pair is named \code{_calculated} to keep it distinct from
#'     \code{topTable}'s \code{CI.L}/\code{CI.R}, which these tables do not
#'     carry: for a \code{dream} fit those bound every contrast by the first
#'     contrast's degrees of freedom.}
#'     \item{degrees_of_freedom}{numeric; the per-feature residual degrees of
#'     freedom used to shrink that feature's variance. Note this is NOT the
#'     degrees of freedom the p-value was computed from, which is the
#'     per-contrast Satterthwaite value plus the prior; p-values cannot be
#'     recomputed from this column.}
#'     \item{logLik}{numeric; log likelihood of differential expression.}
#'     \item{AveExpr}{numeric; mean of all sample-level values for that
#'     feature.}
#'     \item{methylation_diff}{numeric; methylation difference.}
#'     \item{t}{numeric; moderated t-statistic.}
#'     \item{z.std}{numeric; standard normal equivalent of the t-statistic
#'     (z-scores).}
#'     \item{p_value}{numeric; p-value.}
#'     \item{adj_p_value}{numeric; p-values adjusted within each combination of
#'     tissue, assay, platform, and contrast using the Benjamini-Hochberg method
#'     to control the false discovery rate.}
#'   }
#'
#' @author Tyler Sagendorf Christopher Jin
#'
#' @importFrom data.table setorderv setcolorder rbindlist
#' @importFrom utils data
#'
#' @export load_differential_analysis
#'
#' @examples
#' DA_list <- load_differential_analysis() # default behavior
#'
#' # Structure of a single object
#' str(DA_list[["adipose"]][["prot-pr"]])
#'
#' # Un-nest list
#' DA_list <- unlist(DA_list, recursive = FALSE)
#' names(DA_list)
#'
#' # Include epigen data
#' \dontrun{
#' DA_list <- load_differential_analysis(epigen = TRUE)
#' }
#'

load_differential_analysis <- function(selected_omes = "all",
                                       selected_tissues = "all",
                                       single_matrix = FALSE,
                                       epigen = FALSE,
                                       repo_local_dir = NULL,
                                       combine_with_featgene = FALSE,
                                       verbose = TRUE,
                                       load_clinical = FALSE) {
  if (!is.null(repo_local_dir)) {
    message("`repo_local_dir` is no longer needed and is ignored: epigenomics ",
            "differential analysis is read from the public CloudFront release.")
  }

  selected_tissues <- match.arg(
    arg = selected_tissues,
    choices = c("all", "adipose", "blood", "muscle"),
    several.ok = TRUE
  )

  #-----here I basically just make sure that if any metab platform is listed,
  #all metab is loaded, to support differences in platform specific loading
  #
  # metab-t-clinical is exempt. Differential analysis combines the metabolomics
  # platforms into one *_METAB_DA table per tissue, but clinical metabolomics is
  # kept out of it and published separately as BLOOD_METAB_T_CLINICAL_DA. Folding
  # it into "metab" would quietly return the combined table instead of the one
  # that was asked for.
  metab_platforms <- grepl("metab", selected_omes) &
    !selected_omes %in% clinical_ome_list()
  if (any(metab_platforms)) {
    selected_omes = c(selected_omes[!metab_platforms], "metab")
  }

  selected_omes <- match.arg(
    arg = selected_omes,
    choices = c(
      "all", "transcript-rna-seq", "prot-pr", "prot-ph", "prot-ol",
      "metab", "epigen-atac-seq", "epigen-methylcap-seq",
      clinical_ome_list()
    ),
    several.ok = TRUE
  )

  if ("all" %in% selected_tissues) {
    selected_tissues <- c("adipose", "blood", "muscle")
  }

  # Handle epigen omes requested (explicitly or via "all") while epigen is off
  epigen_omes <- c("epigen-atac-seq", "epigen-methylcap-seq")
  requested_epigen <- intersect(selected_omes, epigen_omes)
  if (!epigen & length(requested_epigen) > 0 &
      all(selected_omes %in% epigen_omes)) {
    stop(
      "You've requested only epigenetic omes (",
      paste(requested_epigen, collapse = ", "),
      ") but `epigen = FALSE`. Set `epigen = TRUE` to load epigenetic data."
    )
  }
  if (verbose & !epigen &
      (length(requested_epigen) > 0 | "all" %in% selected_omes)) {
    message(
      "You've requested one or more epigenetic omes (via explicit selection ",
      "or \"all\") but `epigen = FALSE`, so epigenetic data will be skipped. ",
      "Set `epigen = TRUE` to load epigenetic data."
    )
  }

  if ("all" %in% selected_omes) {
    selected_omes <- c(
      "transcript-rna-seq", "prot-pr", "prot-ph", "prot-ol", "metab",
      "epigen-atac-seq", "epigen-methylcap-seq",
      clinical_ome_list()
    )
  }

  # Clinical chemistry is opt-in, the same way epigenomics is. Applied after
  # both expansions so it governs "all" and a named request alike.
  if (!load_clinical) {
    dropped <- base::intersect(selected_omes, clinical_ome_list())
    remaining <- base::setdiff(selected_omes, clinical_ome_list())
    # Asking only for what the gate removes leaves nothing to load, and an empty
    # selection surfaces further down as a data.table error about a missing
    # column. Say what actually happened.
    if (length(dropped) && !length(remaining)) {
      stop("You've requested only clinical omes (",
           paste(dropped, collapse = ", "),
           ") but `load_clinical = FALSE`. Set `load_clinical = TRUE` to load ",
           "clinical chemistry.")
    }
    selected_omes <- remaining
    if (length(dropped) && verbose) {
      message("Clinical omes (", paste(dropped, collapse = ", "),
              ") are skipped; set `load_clinical = TRUE` to include them.")
    }
  }

  if (epigen) {
    selected_omes_epigen <- selected_omes[selected_omes %in%
                                            c("epigen-atac-seq",
                                              "epigen-methylcap-seq")]
    if(verbose){
      message("You've elected to load in the epigenetic data too. These file sizes are significantly larger and will be downloaded from ", .AWS_EPIGEN_DA_URL, " on every call. This loading can be quite slow.")
    }
  }

  if(verbose & "metab" %in% selected_omes){
    message("Please remember that the lowest CV Metabolite is chosen and the
            relevant refmet name is used. If you're not able to find your desired
            metabolite, look through the METABOLOMICS_CVS object for the relevant
            refmet/feature name.")
  }
  # Split epigen platforms from non-epigen, load epigen via old functionality
  selected_omes <- selected_omes[!selected_omes %in%
                                   c("epigen-atac-seq",
                                     "epigen-methylcap-seq")]

  DA_files <- data(package = "MotrpacHumanPreSuspensionAnalysis")
  DA_files <- DA_files[["results"]][, "Item"]
  DA_files <- DA_files[grepl("_DA$", DA_files)]

  tissues <- tolower(sub("\\_.*", "", DA_files))

  omes <- sub("^[^_]+_(.*)_DA$", "\\1", DA_files)
  # gsub, not sub: an ome name carries one underscore per hyphen, so replacing
  # only the first leaves anything with three or more parts malformed and it
  # then matches no request. BLOOD_METAB_T_CLINICAL_DA derived as
  # "metab-t_clinical" rather than "metab-t-clinical" and was unreachable; the
  # epigen tables had the same defect, masked only because they are split off
  # above and downloaded. load_summary_stats() and load_qc() already
  # gsub.
  omes <- gsub("_", "-", tolower(omes))
  omes[omes == "trnscrpt"] <- "transcript-rna-seq"

  new_names <- structure(
    .Data = paste0(tissues, ".", omes),
    names = DA_files
  )

  keep <- (tissues %in% selected_tissues) & (omes %in% selected_omes)
  new_names <- new_names[keep]

  # Load DA results into a list. Only works because of lazy loading
  out <- vector(mode = "list", length = length(new_names))
  names(out) <- as.character(new_names)

  # Objects are returned as they are stored. Every metabolomics table, clinical
  # chemistry included, reads assay = "metab" and names its platform in the
  # `platform` column, so nothing is relabelled on read.
  #
  # A previous version rewrote BLOOD_METAB_T_CLINICAL_DA's assay to
  # "metab-t-clinical", because the summary statistics of the day put the
  # platform in `assay` and the two tiers therefore disagreed. They no longer
  # do: BLOOD_METAB_T_CLINICAL_SUM_STATS carries the same assay/platform pair
  # this object does.
  #
  # What that means for callers: `assay` names the assay family, not the
  # platform. Clinical chemistry and the research platforms both read "metab",
  # and five analytes — Cortisol, Glycerol, KET, NEFA and Glucose — exist on
  # both, so a key that must tell them apart has to include `platform`.
  # (tissue, assay, feature_id) alone selects two rows for those five, silently.
  # `load_clinical = FALSE` is the default, so clinical rows only arrive when
  # they were asked for.
  for (i in seq_along(new_names)) {
    out[[i]] <- eval(parse(text = names(new_names[i])))
  }

  if (epigen) {
    epi_list <- .load_DA_from_AWS(selected_tissues = selected_tissues,
                                  selected_omes = selected_omes_epigen)

    epi_list <- unlist(epi_list, recursive = FALSE)
    epi_list <- .process_raw_DA(epi_list)
    out <- c(out, epi_list)
  }

  if (combine_with_featgene) {
    out <- lapply(out, function(xi) {
      cols <- colnames(xi)

      xi <- merge(
        x = xi,
        y = MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE,
        by = c("assay", "feature_id"),
        all.x = TRUE,
        all.y = FALSE
      )

      setcolorder(xi, neworder = cols) # reset column order
      setcolorder(
        x = xi,
        neworder = setdiff(
          x = colnames(MotrpacHumanPreSuspensionAnalysis::HUMAN_FEATURE_TO_GENE),
          y = cols
        ),
        after = "feature_id"
      )

      # Remove columns with only missing values
      keep_cols <- vapply(xi, function(col_i) any(!is.na(col_i)), logical(1L))

      xi <- xi[, which(keep_cols), with = FALSE]
    })
  }

  if (single_matrix) {
    out <- rbindlist(l = out, use.names = TRUE, fill = TRUE)

    setorderv(x = out, cols = "contrast", order = 1L)
  } else {
    # Nest by tissue
    tissues <- sub("\\..*$", "", names(out))
    names(out) <- sub(".*\\.", "", names(out))
    out <- split(do.call(list, out), tissues)
  }

  return(out)
}


## Internal functions ----------------------------------------------------------

#' @title Process Raw DA Results
#'
#' @description Converts each \code{data.frame} in a list of DA results to a
#'   \code{data.table}, converts character columns to factors, reorders columns,
#'   and sets the key.
#'
#' @param DA_list a named list of DA results (individual \code{data.frame}
#'   objects).
#'
#' @returns A modified version of \code{DA_list} where each list element is a
#'   keyed \code{data.table}. The object will use significantly less memory.
#'
#' @importFrom dplyr %>% left_join select arrange mutate across any_of relocate everything
#' @importFrom data.table as.data.table := setcolorder setorderv setkeyv copy setDT
#'
#' @author Tyler Sagendorf Christopher Jin
#'
#'
#' @noRd

.process_raw_DA <- function(DA_list) {
  for (name_i in names(DA_list)) {
    dt <- copy(DA_list[[name_i]])
    setDT(x = dt)

    # Add contrast information
    dt <- merge(
      x = dt, y = MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER, by = "contrast",
      all.x = TRUE, all.y = FALSE
    )

    contrast_levels <- levels(x = MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER[["contrast"]])

    dt[, `:=`(
      tissue = sub("\\..*$", "", name_i),
      assay = sub("^.*\\.", "", name_i)
    )]

    dt[, `:=`(
      tissue = as.character(tissue),
      assay = as.character(assay),
      contrast = factor(x = contrast, levels = contrast_levels),
      feature_id = as.character(feature_id),
      full_model = as.factor(full_model)
    )]

    dt[, contrast_order := NULL]

    if ("platform" %in% colnames(dt)) {
      dt[, platform := as.factor(platform)]
    }

    dt[, contrast := droplevels(contrast)]

    # Reorder columns
    new_order <- intersect(
      x = c(
        "tissue", "assay", "platform", "full_model",
        colnames(MotrpacHumanPreSuspensionAnalysis::CONTRAST_CONVERTER)
      ),
      y = colnames(dt)
    )
    setcolorder(x = dt, neworder = new_order)

    if ("z.std" %in% colnames(dt)) {
      setcolorder(x = dt, neworder = "z.std", before = "p_value")
    }

    # Reorder rows
    keys <- intersect(
      x = c(
        "full_model", "contrast", "platform",
        "p_value", "feature_id"
      ),
      y = colnames(dt)
    )
    setkeyv(x = dt, cols = keys)

    DA_list[[name_i]] <- dt
  }

  return(DA_list)
}



#' @title Prioritize redundancies using CVs, Uses refmet names
#'
#' @param file_loaded a data frame of metabolomics differential abundance with
#' the column name structure specified above.
#'
#' @description
#' There are also some situations where multiple feature_ids within one platform
#' map to one given refmet name (e.g. metabolite_isomer1::refmet name & metabolite_isomer2::refmet name).
#' For situations like this, I will be also just choosing lowest CV amongst the
#' options and replacing the name with the refmet name.
#'
#' Every time a load_qc or load_da is called, features should be using refmet
#' names and all spelling/naming inconsistencies will be using refmet stuff.
#'
#' @returns a data frame with no redudant metabolites at a refmet level.
#'
#' @noRd
.prioritize_metab_by_cv_da = function(file_loaded,
                                      tissue,
                                      ome){
  lowest_cv_check = MotrpacHumanPreSuspensionAnalysis::METABOLOMICS_CVS %>%
    dplyr::filter(tissue == !!tissue,
                  assay == ome) %>%
    dplyr::select(feature_id, refmet_name, lowest_CV)
  file_loaded = file_loaded %>%
    dplyr::left_join(., lowest_cv_check, by = "feature_id")  %>%
    dplyr::filter(lowest_CV == "yes") %>%
    dplyr::group_by(contrast) %>%
    dplyr::mutate(adj_p_value = p.adjust(p_value, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(feature_id = refmet_name) %>%
    dplyr::select(-c(refmet_name, lowest_CV)) %>%
    dplyr::filter(!is.na(feature_id))
  return(file_loaded)
}

