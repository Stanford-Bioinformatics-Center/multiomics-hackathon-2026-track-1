# Exhaustive screen of all 28 Table 3 entries in available human aggregate data.
# Base R. Output preserves negative and uncertain mappings and all relevant contrasts.
args<-commandArgs(trailingOnly=TRUE)
p<-if(length(args)) args[1] else "."
out<-file.path(p,"table3");dir.create(out,showWarnings=FALSE)
candidates<-read.csv(file.path(out,"table3_candidates.csv"),stringsAsFactors=FALSE,na.strings="")
sha<-"535b4044e7417413de471104c619120337602b77"
objects<-c('BLOOD_TRNSCRPT_DA','BLOOD_PROT_OL_DA','BLOOD_METAB_DA',
 'BLOOD_METAB_T_CLINICAL_DA','BLOOD_PROT_CLINICAL_DA',
 'MUSCLE_TRNSCRPT_DA','MUSCLE_PROT_PR_DA','MUSCLE_PROT_PH_DA','MUSCLE_METAB_DA',
 'ADIPOSE_TRNSCRPT_DA','ADIPOSE_PROT_PR_DA','ADIPOSE_PROT_PH_DA','ADIPOSE_METAB_DA')
readobj<-function(obj){
 f<-file.path(p,paste0(obj,'.rda'))
 if(!file.exists(f)) download.file(paste0('https://raw.githubusercontent.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/',sha,'/data/',obj,'.rda'),f,mode='wb',quiet=TRUE)
 e<-new.env();load(f,e);as.data.frame(e[[obj]])
}
mp<-readobj('HUMAN_FEATURE_TO_GENE')
mm<-unique(mp[,c('assay','feature_id','gene_symbol')])
split_values<-function(x)if(is.na(x)||!nzchar(x))character(0)else strsplit(x,';',fixed=TRUE)[[1]]
matches<-list();inventory<-list();k<-0
for(obj in objects){
 d<-readobj(obj)
 d<-d[d$contrast_type %in% c('exercise_with_controls','Endur_vs_Resist'),]
 if(!'platform'%in%names(d))d$platform<-as.character(d$assay)
 inventory[[obj]]<-data.frame(object=obj,tissue=unique(d$tissue),assay=unique(d$assay),tested_features=length(unique(d$feature_id)),rows=nrow(d))
 for(i in seq_len(nrow(candidates))){
  c<-candidates[i,];s<-d[FALSE,];s$matched_gene<-character(0)
  if(c$mapping_status=='gene'){
   gm<-mm[!is.na(mm$gene_symbol)&mm$gene_symbol %in% split_values(c$genes)&mm$assay %in% unique(d$assay),]
   names(gm)[names(gm)=='gene_symbol']<-'matched_gene'
   if(nrow(gm))gm<-aggregate(matched_gene~assay+feature_id,data=gm,FUN=function(v)paste(sort(unique(v)),collapse=';'))
   s<-merge(d,gm,by=c('assay','feature_id'))
  }else if(any(grepl('metab',d$assay))){
   s<-d[as.character(d$feature_id) %in% split_values(c$metabolite_candidates),]
   s$matched_gene<-rep(NA_character_,nrow(s))
  }
  if(nrow(s)){
   s$candidate_id<-c$candidate_id;s$candidate_name<-c$name;s$mapping_status<-c$mapping_status
   s$source_object<-obj
   s$evidence_layer<-ifelse(s$assay=='transcript-rna-seq','RNA',ifelse(s$assay=='prot-ph','phosphosite',ifelse(grepl('metab',s$assay),'metabolite','protein')))
   s$interpretation_scope<-ifelse(s$tissue=='blood'&s$evidence_layer %in% c('protein','metabolite'),'circulating_abundance','tissue_or_blood_cell_context')
   s$identity_resolved<-s$mapping_status!='ambiguous'
   s$direction<-ifelse(s$logFC>0,'up',ifelse(s$logFC<0,'down','zero'))
   s$review_direction<-c$review_acute_plasma_direction
   # Direction comparisons only apply to circulating abundance, and no change
   # is never inferred from nonsignificance.
   s$direction_comparison<-ifelse(s$interpretation_scope!='circulating_abundance' | !s$identity_resolved,'not_applicable',ifelse(s$adj_p_value>=.05,'no_detected_change',ifelse(s$direction %in% split_values(c$review_acute_plasma_direction),'consistent_with_review','opposite_to_review')))
   cols<-c('candidate_id','candidate_name','mapping_status','identity_resolved','source_object','tissue','assay','platform','evidence_layer','interpretation_scope','feature_id','matched_gene','contrast_type','contrast_category','Timepoint','logFC','CI.L_calculated','CI.R_calculated','p_value','adj_p_value','direction','review_direction','direction_comparison')
   k<-k+1;matches[[k]]<-s[,cols]
  }
 }
 cat(obj,'checked\n')
}
res<-do.call(rbind,matches);rownames(res)<-NULL
# Deduplicate exact feature/gene mapping duplicates; each assay feature remains separate.
res<-unique(res)
write.csv(res,file.path(out,'all_candidate_results.csv'),row.names=FALSE)
write.csv(do.call(rbind,inventory),file.path(out,'assay_inventory.csv'),row.names=FALSE)
cir<-res[res$interpretation_scope=='circulating_abundance'&res$contrast_type=='exercise_with_controls',]
# Conservative family-wide correction over all resolved measured circulating
# candidate-by-assay-by-time-by-mode tests. Retain source FDR as primary metadata.
cir$table3_screen_bh<-NA_real_;valid<-cir$identity_resolved&!is.na(cir$p_value)
cir$table3_screen_bh[valid]<-p.adjust(cir$p_value[valid],method='BH')
write.csv(cir,file.path(out,'plasma_results.csv'),row.names=FALSE)
sm<-list();layers<-c('blood|prot-ol','blood|metab','blood|transcript-rna-seq',
 'muscle|transcript-rna-seq','muscle|prot-pr','muscle|prot-ph','muscle|metab',
 'adipose|transcript-rna-seq','adipose|prot-pr','adipose|prot-ph','adipose|metab')
for(i in seq_len(nrow(candidates))){
 c<-candidates[i,];a<-res[res$candidate_id==c$candidate_id&res$contrast_type=='exercise_with_controls',]
 b<-cir[cir$candidate_id==c$candidate_id,]
 c$plasma_features<-length(unique(paste(b$assay,b$platform,b$feature_id)))
 c$plasma_test_count<-nrow(b)
 c$plasma_status<-if(c$mapping_status=='ambiguous')'identity_unresolved'else if(!nrow(b))'absent_from_analyzed_tables'else if(any(b$adj_p_value<.05))'response_detected'else'present_no_source_FDR_hit'
 c$plasma_source_fdr_hits<-sum(b$adj_p_value<.05)
 c$plasma_screen_bh_hits<-sum(b$table3_screen_bh<.05,na.rm=TRUE)
 c$plasma_hit_directions<-paste(sort(unique(b$direction[b$adj_p_value<.05])),collapse=';')
 c$plasma_min_source_fdr<-if(nrow(b))min(b$adj_p_value)else NA_real_
 c$plasma_min_screen_bh<-if(any(!is.na(b$table3_screen_bh)))min(b$table3_screen_bh,na.rm=TRUE)else NA_real_
 for(l in layers){
  v<-strsplit(l,'|',fixed=TRUE)[[1]];z<-a[a$tissue==v[1]&a$assay==v[2],]
  c[[l]]<-if(!nrow(z))'absent'else if(c$mapping_status=='ambiguous')'identity_unresolved'else if(!any(z$adj_p_value<.05))'tested_NS'else paste(sort(unique(z$direction[z$adj_p_value<.05])),collapse=';')
 }
 sm[[i]]<-c
}
summary<-do.call(rbind,sm)
stopifnot(nrow(summary)==28,all(res$adj_p_value>=0&res$adj_p_value<=1))
write.csv(summary,file.path(out,'coverage_summary.csv'),row.names=FALSE)
print(summary[,c('name','plasma_status','plasma_features','plasma_source_fdr_hits','plasma_screen_bh_hits','plasma_hit_directions')],row.names=FALSE)
cat('Candidate result rows:',nrow(res),'\nCirculating hypothesis-family tests:',sum(valid),'\n')
