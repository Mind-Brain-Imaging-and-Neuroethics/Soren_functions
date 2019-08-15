mxmnrfit <- function(filein,frmla,fileout){
  require('RJSONIO')
  require('brms')
  
  library('brms')
  library('RJSONIO')
  
  mydata <- read.csv(filein)
  
  #mdl <- MCMCglmm(fixed_frmla,random=random_frmla,data=mydata)
  
  p <- set_prior('cauchy(0,2.5)',class='b')
  
  mdl <- brm(frmla,data=mydata,family=categorical(),prior=p,iter=100000,control=list(max_treedepth=15))
  
  mdl_json <- toJSON(mdl)
  
  smry <- summary(mdl)
  
  smry_json <- toJSON(smry)
  
  write(mdl_json,file=paste(fileout,'_model.json',sep=''))
  
  write(smry_json,file=paste(fileout,'_summary.json',sep=''))
  
}
  
  
  