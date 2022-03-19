isolate <- function(x, node = NULL, depth = 0){
  Qx <- quantile(x, probs = seq(0, 1, 0.01))
  if(is.null(node)){
    node = sample(x, 1)
  }else{
    node = sample
  }

  while(depthMax = F){
    ZeroNode = sample(x, 1)
    Zstat = min(which(Qx >= ZeroNode))-1
    MaxNode = max(x) - ZeroNode
  }
}

#bin by quantile
#sample random value from data
#subset the data to only values which are > sample
#keep doing this until you have no values left
#keep track of 'H' many time you did this
#repeat for 'n' number of times
#values which have big H are more likely outliers