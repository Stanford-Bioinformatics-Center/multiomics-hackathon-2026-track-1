# Human c2.0 discovery screen. Recompute BH before selecting genes.
# RNA and total protein support are kept separate; phosphosites are not abundance evidence.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2)
input <- normalizePath(args[1], mustWork=TRUE)
output <- args[2]; dir.create(output, recursive=TRUE, showWarnings=FALSE)
read_object <- function(name) {
  env <- new.env(parent=emptyenv()); load(file.path(input, paste0(name, '.rda')), env)
  as.data.frame(env[[name]])
}
mapping <- read_object('HUMAN_FEATURE_TO_GENE')
mapping <- unique(mapping[!is.na(mapping$gene_symbol) & mapping$assay %in% c('prot-ol','prot-pr','transcript-rna-seq'),
                          c('assay','feature_id','gene_symbol','uniprot')])
collapse <- function(x) { x <- as.character(x); paste(sort(unique(x[!is.na(x) & nzchar(x)])), collapse=';') }
gene_map <- aggregate(gene_symbol ~ assay + feature_id, mapping, collapse)
protein_map <- aggregate(uniprot ~ assay + feature_id, mapping, collapse)
gene_map <- merge(gene_map, protein_map, by=c('assay','feature_id'), all.x=TRUE)
gene_map$ambiguous_gene <- grepl(';', gene_map$gene_symbol, fixed=TRUE)
plasma_genes <- unique(gene_map$gene_symbol[gene_map$assay == 'prot-ol' & !gene_map$ambiguous_gene])
objects <- c('BLOOD_PROT_OL_DA','BLOOD_TRNSCRPT_DA','MUSCLE_TRNSCRPT_DA','ADIPOSE_TRNSCRPT_DA',
             'MUSCLE_PROT_PR_DA','ADIPOSE_PROT_PR_DA')
tables <- list(); audits <- list(); inventories <- list(); models <- list()
for (object in objects) {
  dat <- read_object(object)
  if (!'platform' %in% names(dat)) dat$platform <- as.character(dat$assay)
  inventories[[object]] <- data.frame(source_object=object, tissue=unique(dat$tissue), assay=unique(dat$assay),
                                    all_features=length(unique(dat$feature_id)), all_contrast_rows=nrow(dat))
  dat <- dat[dat$contrast_type %in% c('exercise_with_controls','Endur_vs_Resist'), ]
  stopifnot(all(is.finite(dat$p_value)), all(dat$p_value >= 0 & dat$p_value <= 1),
            !anyDuplicated(dat[,c('tissue','assay','platform','contrast','feature_id')]))
  family <- interaction(dat$tissue, dat$assay, dat$platform, dat$contrast, drop=TRUE)
  dat$bh_recomputed <- ave(dat$p_value, family, FUN=function(x) p.adjust(x, method='BH'))
  delta <- abs(dat$bh_recomputed - dat$adj_p_value)
  stopifnot(all(is.finite(delta)), max(delta) < 1e-10)
  dat$n_features_in_assay <- ave(dat$p_value, family, FUN=length)
  for (ix in split(seq_len(nrow(dat)), family)) {
    q <- dat[ix, ]; first <- q[1, ]
    audits[[length(audits)+1]] <- data.frame(source_object=object, tissue=first$tissue, assay=first$assay,
      platform=first$platform, contrast=first$contrast, contrast_category=first$contrast_category,
      Timepoint=first$Timepoint, n_features=nrow(q), source_FDR_hits=sum(q$adj_p_value < .05),
      source_FDR_up=sum(q$adj_p_value < .05 & q$logFC > 0), max_BH_difference=max(abs(q$bh_recomputed-q$adj_p_value)))
  }
  models[[object]] <- data.frame(source_object=object, full_model=unique(dat$full_model))
  dat$source_object <- object
  keep <- c('source_object','tissue','assay','platform','feature_id','contrast','contrast_type','contrast_category',
            'Timepoint','logFC','CI.L_calculated','CI.R_calculated','p_value','adj_p_value','bh_recomputed','n_features_in_assay')
  dat <- merge(dat[,keep], gene_map, by=c('assay','feature_id'), all.x=TRUE)
  if (object == 'BLOOD_PROT_OL_DA') {
    # Keep every Olink feature, including unmapped ones, for Python BH cross-checks.
    write.csv(dat, file.path(output,'plasma_full_background.csv'), row.names=FALSE)
  } else {
    dat <- dat[!is.na(dat$gene_symbol) & dat$gene_symbol %in% plasma_genes & !dat$ambiguous_gene, ]
  }
  tables[[object]] <- dat
  cat(object, ': retained', nrow(dat), 'gene-matched rows after full-assay BH validation\n')
}
con <- gzfile(file.path(output,'mapped_source_results.csv.gz'),'wt')
write.csv(do.call(rbind,tables), con, row.names=FALSE); close(con)
write.csv(do.call(rbind,audits),file.path(output,'bh_audit.csv'),row.names=FALSE)
write.csv(do.call(rbind,inventories),file.path(output,'assay_inventory.csv'),row.names=FALSE)
write.csv(do.call(rbind,models),file.path(output,'source_models.csv'),row.names=FALSE)
write.csv(gene_map[gene_map$assay == 'prot-ol',],file.path(output,'plasma_gene_mapping.csv'),row.names=FALSE)
