mxmnrfit <- function(filein,frmla,fileout){
  require('RJSONIO')
  require('brms')
  
  library('brms')
  library('RJSONIO')
  
  mydata <- read.csv(filein)
  
  #mdl <- MCMCglmm(fixed_frmla,random=random_frmla,data=mydata)
  
  mdl <- brm(frmla,data=mydata,family=categorical())
  
  options(expressions=10000)
  
  mdl_json <- toJSON(mdl)
  
  smry <- summary(mdl)
  
  smry_json <- toJSON(smry)
  
  write(mdl_json,file=paste(fileout,'_model.json',sep=''))
  
  write(smry_json,file=paste(fileout,'_summary.json',sep=''))
  
}
  
  
  