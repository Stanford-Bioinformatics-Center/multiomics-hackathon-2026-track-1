# Published human c2.0 summaries; no individual-level refitting.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
input <- normalizePath(args[1],mustWork=TRUE); output <- args[2]
dir.create(output,recursive=TRUE,showWarnings=FALSE)
read_object <- function(n) {
  e <- new.env(parent=emptyenv()); load(file.path(input,paste0(n,'.rda')),e)
  as.data.frame(e[[n]])
}
mapping <- read_object('HUMAN_FEATURE_TO_GENE')
map <- unique(mapping[!is.na(mapping$gene_symbol) & mapping$gene_symbol %in% c('WARS1','WARS','WRS','IFI53'),
                      c('assay','feature_id','gene_symbol','entrez_gene','ensembl_gene','uniprot')])
objects <- c('BLOOD_PROT_OL_DA','BLOOD_TRNSCRPT_DA','MUSCLE_TRNSCRPT_DA','ADIPOSE_TRNSCRPT_DA',
             'MUSCLE_PROT_PR_DA','ADIPOSE_PROT_PR_DA','MUSCLE_PROT_PH_DA','ADIPOSE_PROT_PH_DA')
results <- list(); audits <- list(); coverage <- list()
for (object in objects) {
  d <- read_object(object)
  if(!'platform' %in% names(d)) d$platform <- as.character(d$assay)
  lookup <- unique(map[map$assay %in% unique(d$assay),c('assay','feature_id','gene_symbol')])
  stopifnot(!anyDuplicated(lookup[,c('assay','feature_id')]))
  hits <- d$feature_id %in% lookup$feature_id
  coverage[[object]] <- data.frame(source_object=object,tissue=unique(d$tissue),assay=unique(d$assay),
    source_features=length(unique(d$feature_id)),mapped_features=length(unique(lookup$feature_id)),
    analyzed_WARS1_features=length(unique(d$feature_id[hits])),source_WARS1_rows=sum(hits))
  if(!any(hits)) {cat(object,': no analyzed WARS1 features\n');next}
  # Full feature universe, including all within-arm and control comparisons.
  family <- interaction(d$tissue,d$assay,d$platform,d$contrast,drop=TRUE)
  d$bh_recomputed <- ave(d$p_value,family,FUN=function(p)p.adjust(p,'BH'))
  stopifnot(all(is.finite(d$p_value)),max(abs(d$bh_recomputed-d$adj_p_value))<1e-10)
  for (ix in split(seq_len(nrow(d)),family)) {
    z <- d[ix,];a <- z[1,]
    audits[[length(audits)+1]] <- data.frame(source_object=object,tissue=a$tissue,assay=a$assay,
      platform=a$platform,contrast=a$contrast,contrast_type=a$contrast_type,contrast_category=a$contrast_category,
      Timepoint=a$Timepoint,n_features=nrow(z),max_BH_difference=max(abs(z$bh_recomputed-z$adj_p_value)))
  }
  h <- merge(d[hits,],lookup,by=c('assay','feature_id'))
  h$source_object <- object
  cols <- c('source_object','tissue','assay','platform','feature_id','gene_symbol','contrast','contrast_short',
            'contrast_type','contrast_category','Timepoint','full_model','logFC','CI.L_calculated','CI.R_calculated',
            'p_value','adj_p_value','bh_recomputed')
  results[[object]] <- h[,cols]
  cat(object,':',nrow(h),'WARS1 source contrasts\n')
}
write.csv(do.call(rbind,results),file.path(output,'wars1_all_source_contrasts.csv'),row.names=FALSE)
write.csv(do.call(rbind,audits),file.path(output,'bh_audit.csv'),row.names=FALSE)
write.csv(do.call(rbind,coverage),file.path(output,'wars1_assay_coverage.csv'),row.names=FALSE)
write.csv(map,file.path(output,'wars1_feature_mapping.csv'),row.names=FALSE)

# Two literature-motivated metabolites for post hoc mechanistic context, not a pathway test.
# Export exactly Tryptophan and Kynurenine; do not select isotope internal standards.
metabolites <- c('tryptophan','kynurenine')
met_results <- list(); met_checks <- list()
for (object in c('BLOOD_METAB_DA','MUSCLE_METAB_DA','ADIPOSE_METAB_DA')) {
  d <- read_object(object)
  d <- d[d$contrast_type %in% c('exercise_with_controls','Endur_vs_Resist'),]
  h <- d[tolower(trimws(as.character(d$feature_id))) %in% metabolites,]
  h$source_object <- object
  met_results[[object]] <- h
  # Do not force-match the source BH family: the available metabolite table does
  # not reproduce original adjusted p-values with either of these groupings.
  for (fields in list(c('tissue','assay','platform','contrast'),c('tissue','assay','contrast'))) {
    family <- do.call(interaction,c(d[fields],list(drop=TRUE)))
    q <- ave(d$p_value,family,FUN=function(p)p.adjust(p,'BH'))
    met_checks[[length(met_checks)+1]] <- data.frame(source_object=object,
      attempted_family=paste(fields,collapse=';'),max_difference=max(abs(q-d$adj_p_value)),
      interpretation='Source q retained; this grouping is not an independent reproduction of source correction')
  }
}
write.csv(do.call(rbind,met_results),file.path(output,'tryptophan_kynurenine_source_context.csv'),row.names=FALSE)
write.csv(do.call(rbind,met_checks),file.path(output,'metabolite_BH_reproduction_limits.csv'),row.names=FALSE)
cat('Saved separate exploratory metabolite context, retaining original source q values.\n')
