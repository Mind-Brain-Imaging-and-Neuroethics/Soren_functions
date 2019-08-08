mxmnr_compare <- function(filein,frmla1,frmla2,fileout){
  require('RJSONIO')
  require('brms')
  
  library('brms')
  library('RJSONIO')
  
  add_criterion <- brms::add_criterion
  
  mydata <- read.csv(filein)
  
  p <- set_prior('cauchy(0,2.5)',class='b')
  
  mdl1 <- brm(frmla1,data=mydata,family=categorical(),iter=2000,prior=p,control=list(max_treedepth=15,adapt_delta=0.95))
  mdl1 <- add_criterion(mdl1,criterion='kfold')
  
  tmp <- summary(mdl1)
  smry1 <- list(fixed=tmp$fixed,random=tmp$random,samples=mdl1$fit@sim$samples,kfold=mdl1$kfold)
  smry1_json <- toJSON(smry1)
  
  mdl2 <- brm(frmla2,data=mydata,family=categorical(),iter=2000,prior=p,control=list(max_treedepth=15,adapt_delta=0.95))
  mdl2 <- add_criterion(mdl2,criterion='kfold')
  
  tmp <- summary(mdl2)
  smry2 <- list(fixed=tmp$fixed,random=tmp$random,samples=mdl2$fit@sim$samples,kfold=mdl2$kfold)
  smry2_json <- toJSON(smry2)
  
  cmpare <- loo_compare(mdl1,mdl2,criterion='kfold')
  cmpare_json <- toJSON(cmpare)
  
  
  write(smry1_json,file=paste(fileout,'_summary1.json',sep=''))
  
  write(smry2_json,file=paste(fileout,'_summary2.json',sep=''))
  
  write(cmpare_json,file=paste(fileout,'_cmpare.json',sep=''))
  
  
  
  
}
  
  
  