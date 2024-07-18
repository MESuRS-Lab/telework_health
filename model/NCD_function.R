
#daily risk 

NCD_function = function(DRF, theoretical = F, bound = 0.7) {
  
  if(!theoretical){
    if(DRF == "LS"){
      
      NCD_rate = approxfun(c(0, 0.25/5, 1/5, 5/5),
                           c((1/(1+exp(-(-3.33)))) / (1/(1+exp(-(-3.33)))),
                             (1/(1+exp(-(-3.33-0.571)))) / (1/(1+exp(-(-3.33)))),
                             (1/(1+exp(-(-3.33-0.997)))) / (1/(1+exp(-(-3.33)))),
                             (1/(1+exp(-(-3.33-1.233)))) / (1/(1+exp(-(-3.33))))))
      
    } else if(DRF == "US"){
      NCD_rate = approxfun(c(0,0.5,1), c(2.33/2.33, 1.71/2.33, 2.08/2.33))
    } else if(DRF == "IU"){
      NCD_rate = approxfun(c(0,0.2,0.5,1), c(49.9/49.9, 52.8/49.9, 53.2/49.9, 47.5/49.9))
    } else stop("Supported non-theoretical DRF are: LS, US, IU")
  } else {
    if(DRF == "LS"){
      NCD_rate = function(val, fbound=bound){(1-fbound)*exp(-10*val)+fbound}
    } else if(DRF == "US"){
      if(bound>1) stop("A bound greater than 1 will give you an IU not an US shape!")
      NCD_rate = function(val, fbound=bound){4*(1-fbound)*(val-0.5)^2+fbound}
    } else if(DRF == "IU"){
      if(bound<1) stop("A bound lower than 1 will give you a US not an IU shape!")
      NCD_rate = function(val, fbound=bound){4*(1-fbound)*(val-0.5)^2+fbound}
    } else if(DRF == "LD"){
      if(bound>1) stop("A bound greater than 1 will give you an LD not an LI shape!")
      NCD_rate = function(val, fbound=bound){(fbound-1)*val+1}
    } else if(DRF == "LI"){
      if(bound<1) stop("A bound lower than 1 will give you an LD not an LI shape!")
      NCD_rate = function(val, fbound=bound){(fbound-1)*val+1}
    } else stop("Supported non-theoretical DRF are: LS, US, IU, LD, LI")
    
  }
  
  return(NCD_rate)
  
}
