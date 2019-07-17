cvtest_matlab <- function(filein,fileout){

require('cvequality')
require('RJSONIO')
  
  library(cvequality)
  
filedata <- read.csv(filein)

testoutput <- with(filedata,asymptotic_test(X,G))

testjson <- toJSON(testoutput)

write(testjson,file=fileout,sep="")

}