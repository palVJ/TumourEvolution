AlleleFrequency = function(mutationsInEachcell){
  
  print(Sys.time())
  NrOfCells = length(mutationsInEachcell)
  AllMutations = unlist(mutationsInEachcell)
  UniqueMutations = unique(AllMutations)
  UniqueMutations = UniqueMutations[2:length(UniqueMutations)]
  count = 0
  
  library(foreach)
  library(doParallel)
  cl = makeCluster(2)
  registerDoParallel(cl)

  vaf = foreach(i = UniqueMutations, .combine = 'c')%dopar%{
   count = length(which(AllMutations == i))
   if(count > 0 && count/NrOfCells > 0.01 ){
     count/NrOfCells
  }
  }


  stopCluster(cl)
 #  max = 0
 #  for(i in 1:length(mutationsInEachcell)){
 # 
 #    temp = mutationsInEachcell[[i]][length(mutationsInEachcell[[i]])]
 #    if(temp > max){
 #      max = temp
 #    }
 #  }
 # 
 #    vaf = c()
 #   for(j in 1:max) {
 #    count = 0
 #   for(k in 1:length(mutationsInEachcell)){
 #     if(length(which(mutationsInEachcell[[k]]==j))> 0){
 #       count = count +1}
 #   }
 #   if(count > 0 && count/length(mutationsInEachcell) > 0.01 ){
 #    vaf[length(vaf)+1] = count/length(mutationsInEachcell)
 #   }
 # }
  
  
  hist(vaf, freq = 0, breaks = 50)
  print(Sys.time())
  return(vaf)
  
}