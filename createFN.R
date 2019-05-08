createFN <- function(ple,siglength){
  require(fArma)
  if (ple < 1) {
    sig = fbmSim(n = siglength, H = (ple+1)/2, method = 'circ', doplot = FALSE,fgn = TRUE)
  } else {
    sig = fbmSim(n = siglength, H = (ple-1)/2, method = 'circ',doplot = FALSE)
    }
    
  write.csv(sig,file = "/Users/Soren/output.csv")
}