options(warn=1,survey.lonely.psu='adjust')
suppressPackageStartupMessages(library(survey))
args=commandArgs(trailingOnly=TRUE);w=args[1];o=file.path(w,'analysis');dir.create(o,showWarnings=FALSE)
d=read.csv(file.path(w,'remote_sources/revision_20260916_inputs/nhanes_2015_2018_raw_merge.csv'))
# IBM XPT numeric zero (eight zero bytes) can be decoded by the source
# pandas XPORT reader as 2^-260 rather than 0. Restore this single sentinel
# in every raw numeric variable; the 1e-12 relative tolerance accommodates
# decimal CSV round trips. Do not use a broad small-value cutoff, and do
# not convert these real zeros (including INDFMPIR=0) to missing values.
# Official GLU_I/J and DEMO_I/J XPT byte-level verification is recorded in
# analysis/nhanes_xpt_zero_verification.json by verify_nhanes_xpt_zeros.py.
xpt_zero_sentinel=2^(-260)
xpt_zero_audit=list()
for(column in names(d)) if(is.numeric(d[[column]])) {
 original=d[[column]]
 zero_artifact=is.finite(original) & abs(original/xpt_zero_sentinel-1)<=1e-12
 d[[column]][zero_artifact]=0
 xpt_zero_audit[[column]]=data.frame(variable=column,
   n_restored_zero=sum(zero_artifact),n_missing_before=sum(is.na(original)),
   n_missing_after=sum(is.na(d[[column]])),n_zero_after=sum(d[[column]]==0,na.rm=TRUE))
}
write.csv(do.call(rbind,xpt_zero_audit),file.path(o,'nhanes_xpt_zero_repair.csv'),row.names=FALSE)
stopifnot(!any(vapply(d[vapply(d,is.numeric,logical(1))],
  function(x)any(is.finite(x)&abs(x/xpt_zero_sentinel-1)<=1e-12),logical(1))))
# The source codebooks' smallest positive fasting weights are >1,000;
# this guard catches unresolved zero-decoding artifacts before selection.
stopifnot(all(is.na(d$WTSAF2YR)|d$WTSAF2YR==0|d$WTSAF2YR>1000))
d$age=d$RIDAGEYR;d$sex=factor(d$RIAGENDR);d$race=factor(d$RIDRETH3);d$pir=d$INDFMPIR
d$smoking=ifelse(d$SMQ020==2,'Never',ifelse(d$SMQ020==1 & d$SMQ040==3,'Former',ifelse(d$SMQ020==1 & d$SMQ040 %in% c(1,2),'Current',NA)))
d$smoking=factor(d$smoking,levels=c('Never','Former','Current'))
d$diabetes=ifelse(d$DIQ010==1,1,ifelse(d$DIQ010 %in% c(2,3),0,NA))
d$bmi=d$BMXBMI;d$tyg=log(d$LBXTR*d$LBXGLU/2);d$wbc=d$LBXWBCSI;d$nlr=d$LBXNEPCT/d$LBXLYPCT;d$crp_unadjusted=d$LBXHSCRP
d$crp=d$crp_unadjusted
ix=d$cycle=='2015-2016' & !is.na(d$crp) & d$crp<=23
d$crp[ix]=0.8695*d$crp[ix]+0.2954
hit=rowSums(cbind(d$MCQ230A==23,d$MCQ230B==23,d$MCQ230C==23),na.rm=TRUE)>0
known_type=rowSums(cbind(d$MCQ230A %in% 10:39,d$MCQ230B %in% 10:39,d$MCQ230C %in% 10:39))>0
d$lung_cancer=ifelse(hit,1,ifelse(d$MCQ220==2 | (d$MCQ220==1 & known_type),0,NA))
d$weight=d$WTSAF2YR/2
# Pregnancy coding is inapplicable in most men and older adults: retain missing.
d$eligible=d$age>=40 & !is.na(d$smoking) & complete.cases(d[,c('tyg','bmi','crp','wbc','nlr')]) & (is.na(d$RIDEXPRG)|d$RIDEXPRG!=1)
d$eligible[is.na(d$eligible)]=FALSE
logvars=c('crp','wbc','nlr');for(v in logvars)d[[paste0('log_',v)]]=log1p(d[[v]])
valid=d$eligible & !is.na(d$weight) & d$weight>0
zfun=function(x){m=weighted.mean(x[valid],d$weight[valid],na.rm=TRUE);s=sqrt(weighted.mean((x[valid]-m)^2,d$weight[valid],na.rm=TRUE));(x-m)/s}
for(v in c('tyg','bmi','log_crp','log_wbc','log_nlr'))d[[paste0('z_',v)]]=zfun(d[[v]])
components=c('z_tyg','z_bmi','z_log_crp','z_log_wbc','z_log_nlr')
d$MIDS_full=rowMeans(d[,components]);d$MIDS_no_crp=rowMeans(d[,setdiff(components,'z_log_crp')]);d$MIDS_no_wbc=rowMeans(d[,setdiff(components,'z_log_wbc')]);d$MIDS_no_nlr=rowMeans(d[,setdiff(components,'z_log_nlr')]);d$metabolic_only=rowMeans(d[,c('z_tyg','z_bmi')])
for(v in c('MIDS_full','MIDS_no_crp','MIDS_no_wbc','MIDS_no_nlr','metabolic_only'))d[[v]]=zfun(d[[v]])
selection=data.frame(step=c('All interviewed','Age 40 or older','Positive fasting weight','Complete five components and smoking','Complete primary covariates'),n=c(nrow(d),sum(d$age>=40,na.rm=TRUE),sum(d$age>=40 & d$weight>0,na.rm=TRUE),sum(valid),sum(valid & complete.cases(d[,c('age','sex','race','pir','smoking','diabetes')]))))
write.csv(selection,file.path(o,'nhanes_selection.csv'),row.names=FALSE)
d=d[!is.na(d$weight)&d$weight>0 & !is.na(d$SDMVSTRA)&!is.na(d$SDMVPSU),]
des=svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~weight,nest=TRUE,data=d)
adult=subset(des,eligible)
cut=as.numeric(coef(svyquantile(~MIDS_full,adult,quantiles=c(.25,.5,.75),ci=FALSE,na.rm=TRUE)))
d$quartile=cut(d$MIDS_full,breaks=c(-Inf,cut,Inf),labels=paste0('Q',1:4),include.lowest=TRUE)
write.csv(d[d$eligible,],file.path(o,'nhanes_analysis_data.csv'),row.names=FALSE)
des=update(des,quartile=d$quartile);adult=subset(des,eligible)
rows=list();k=0
model=function(y,x,scope='Primary',label='Adjusted'){
 de=switch(scope,Primary=adult,Ever_smokers=subset(adult,smoking!='Never'),Age60=subset(adult,age>=60),CRP10=subset(adult,crp<=10),Cycle2015=subset(adult,cycle=='2015-2016'),Cycle2017=subset(adult,cycle=='2017-2018'))
 vars=c(y,x,'age','sex','race','pir','smoking','diabetes','cycle');dd=de$variables;keep=complete.cases(dd[,vars]);de=de[keep,];dd=de$variables
 covs=if(label=='Age_sex')c('age','sex','cycle')else c('age','sex','race','pir','smoking','diabetes','cycle')
 covs=covs[vapply(dd[,covs,drop=FALSE],function(v)length(unique(v))>1,logical(1))]
 formula=as.formula(paste(y,'~',paste(c(x,covs),collapse='+')))
 binary=y=='lung_cancer';n=nrow(dd);events=if(binary)sum(dd[[y]])else NA
 if(binary)return(data.frame(outcome=y,exposure=x,scope=scope,model=label,n=n,events=events,estimate=NA,ci_low=NA,ci_high=NA,p=NA,status='Descriptive only: sparse cases and separated covariate strata'))
 if(n<50)return(data.frame(outcome=y,exposure=x,scope=scope,model=label,n=n,events=events,estimate=NA,ci_low=NA,ci_high=NA,p=NA,status='Insufficient observations'))
 tryCatch({f=svyglm(formula,de,family=if(binary)quasibinomial()else gaussian());s=summary(f)$coefficients;ci=confint(f,parm=x);b=coef(f)[x]; if(binary){b=exp(b);ci=exp(ci)};data.frame(outcome=y,exposure=x,scope=scope,model=label,n=n,events=events,estimate=b,ci_low=ci[1],ci_high=ci[2],p=s[x,4],status='Estimated')},error=function(e)data.frame(outcome=y,exposure=x,scope=scope,model=label,n=n,events=events,estimate=NA,ci_low=NA,ci_high=NA,p=NA,status=conditionMessage(e)))
}
for(scope in c('Primary','Ever_smokers','Age60','CRP10','Cycle2015','Cycle2017'))for(y in c('log_crp','log_wbc','log_nlr','lung_cancer')){
 xs=if(y=='lung_cancer')c('MIDS_full','metabolic_only')else c(paste0('MIDS_no_',sub('log_','',y)),'metabolic_only')
 for(x in xs){k=k+1;rows[[k]]=model(y,x,scope)}
}
r=do.call(rbind,rows);r$q=NA
for(sc in unique(r$scope)){ix=which(r$scope==sc & !is.na(r$p));r$q[ix]=p.adjust(r$p[ix],method='BH')}
write.csv(r,file.path(o,'nhanes_models.csv'),row.names=FALSE)
qrows=list()
for(q in levels(d$quartile)){
 de=subset(adult,quartile==q); dd=de$variables
 rr=data.frame(quartile=q,n=nrow(dd),lung_cases=sum(dd$lung_cancer==1,na.rm=TRUE))
 for(v in c('MIDS_full','metabolic_only','tyg','bmi','crp','wbc','nlr','age'))rr[[v]]=as.numeric(coef(svyquantile(as.formula(paste0('~',v)),de,quantiles=.5,ci=FALSE,na.rm=TRUE)))
 qrows[[q]]=rr
}
write.csv(do.call(rbind,qrows),file.path(o,'nhanes_quartiles.csv'),row.names=FALSE)
candidate=d[d$age>=40,]
write.csv(data.frame(variable=names(candidate),denominator=nrow(candidate),n_missing=vapply(candidate,function(x)sum(is.na(x)),numeric(1))),file.path(o,'nhanes_missingness.csv'),row.names=FALSE)
writeLines(c('NHANES 2015-2018 analysis after IBM XPT zero restoration',
 'Rule: finite raw numeric values matching 2^-260 within relative tolerance 1e-12 are restored to numeric 0; missing values remain missing.',
 'Only strictly positive fasting subsample weights enter the survey design.',
 capture.output(selection),
 'Primary and restricted model denominators:',
 capture.output(unique(r[,c('scope','n','events','status')])),
 'Restored raw numeric zeros:',
 capture.output(do.call(rbind,xpt_zero_audit))),file.path(o,'nhanes_analysis_summary.txt'))
writeLines(capture.output(sessionInfo()),file.path(o,'R_session_info.txt'))
print(selection);print(r[r$scope=='Primary',]);print(do.call(rbind,qrows))
