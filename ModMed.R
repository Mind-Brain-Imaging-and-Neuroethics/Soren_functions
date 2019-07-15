ModMed <- function(filein,fileout){
  require('psych')
  require('mediation')
  require('RJSONIO')
  
  MyData <- read.csv(filein)

  MyData$X.c <- scale(MyData$X,center=TRUE,scale=FALSE)[,]
  MyData$W.c <- scale(MyData$W,center=TRUE,scale=FALSE)[,]

  mediate <- mediation::mediate

  test.modmed <- mediation::test.modmed

  ModMedModel1 <- lm(M ~ X.c*W.c,data=MyData)
  
  ModMedModel2 <- lm(Y ~ X.c*W.c + M,data=MyData)
  
  #a <- toJSON(ModMedModel2)
  
  #summary(Mod.Med.Model.2)

  low.W <- mean(MyData$W)-sd(MyData$W)

  Mod.Med.LowW <- mediate(ModMedModel1,ModMedModel2,covariates=list(W.c=low.W),boot=TRUE,
                        boot.ci.type="bca",sims=2000,treat="X.c",mediator="M")

  high.W <- mean(MyData$W)+sd(MyData$W)

  Mod.Med.highW <- mediate(ModMedModel1,ModMedModel2,covariates=list(W.c=high.W),boot=TRUE,
                        boot.ci.type="bca",sims=2000,treat="X.c",mediator="M")

  Mod.Med.TestW <- mediate(model.m=ModMedModel1, model.y=ModMedModel2, boot = TRUE,  
                            boot.ci.type = "bca", sims = 2000, treat="X.c", mediator="M")

  Mod.Med.Test <- test.modmed(Mod.Med.TestW, covariates.1 = list(W.c = low.W),   
            covariates.2 = list(W.c = high.W), sims = 2000)

  MedModel <- toJSON(Mod.Med.TestW)
  
  ModTest <- toJSON(Mod.Med.Test)
  
  write(MedModel,file=paste(fileout,'_model.json',sep=""))
  
  write(ModTest,file=paste(fileout,'_modtest.json',sep=""))

}


