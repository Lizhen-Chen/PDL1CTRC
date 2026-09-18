library(survival)
d<-read.csv('analysis/sc_clinical_followup.csv',check.names=FALSE)
rows<-list();n<-1
for (f in c('score_sd','non_MPR')) {
 for (label in c('unadjusted','clinical','response_adjusted')) {
  if(f=='non_MPR' && label=='response_adjusted') next
  terms<-c(f,if(label=='clinical')c('age','male','LUSC','stage3') else if(label=='response_adjusted') c('age','male','LUSC','stage3','non_MPR'))
  z<-d[complete.cases(d[,c(terms,'RFS_months','RFS_status')]),]
  fit<-coxph(as.formula(paste('Surv(RFS_months,RFS_status)~',paste(terms,collapse='+'))),data=z,ties='efron',x=TRUE)
  ph<-cox.zph(fit,transform='km')$table
  for(term in rownames(ph)) {rows[[n]]<-data.frame(feature=f,model=label,term=term,n=nrow(z),events=sum(z$RFS_status),chisq=ph[term,1],df=ph[term,2],p=ph[term,3]);n<-n+1}
 }
}
write.csv(do.call(rbind,rows),'analysis/sc_rfs_ph_assumption.csv',row.names=FALSE)
print(do.call(rbind,rows))
