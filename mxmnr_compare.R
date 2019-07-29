mxmnr_compare <- function(filein,frmla1,frmla2,fileout){
  require('RJSONIO')
  require('brms')
  
  library('brms')
  library('RJSONIO')
  
  mydata <- read.csv(filein)
  
  mdl1 <- brm(frmla1,data=mydata,family=categorical())
  
  tmp <- summary(mdl1)
  smry1 <- list(fixed=tmp$fixed,random=tmp$random)
  smry1_json <- toJSON(smry1)
  
  mdl2 <- brm(frmla2,data=mydata,family=categorical())
  
  tmp <- summary(mdl2)
  smry2 <- list(fixed=tmp$fixed,random=tmp$random)
  smry2_json <- toJSON(smry2)
  
  #add_criterion(mdl1,c('waic','loo'))
  #add_criterion(mdl2,c('waic','loo'))
  
  #cmpare <- loo_compare(mdl1,mdl2,criterion=c('waic','loo'))
  #cmpare_json <- toJSON(cmpare)
  
  #smry <- summary(mdl)
  
  #smry_json <- toJSON(smry)
  
  #write(mdl_json,file=paste(fileout,'_model.json',sep=''))
  
  write(smry1_json,file=paste(fileout,'_summary1.json',sep=''))
  
  write(smry2_json,file=paste(fileout,'_summary2.json',sep=''))
  
  #write(cmpare_json,file=paste(fileout,'_cmpare.json',sep=''))
  
  
  
  
}
  
  
  