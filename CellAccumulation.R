BirthAndDeath = function(beta,lambda,n, ploidy = 2){
#Begins with one cancer cell  
#Observe for a total of n generations  
NumberOfCells = 1
MutationsInEachCell = 0
m = 1 #labeled mutation 
for(i in 1:n){
 
  CellsThatDivided = rbinom(1, NumberOfCells, beta)
  #choose which cells divided:
  WhichCellsDivided = sample(1:NumberOfCells,CellsThatDivided,replace = FALSE)
 
  UpdateNumberOfCells = 2*CellsThatDivided

  if(UpdateNumberOfCells == 0){
    return(UpdateNumberOfCells)
  }
  
  #preallocate vector
  UpdateMutationsInEachCell = vector(mode = "list", length = UpdateNumberOfCells)
  for(r in 1:UpdateNumberOfCells){
    UpdateMutationsInEachCell[[r]] = 0
  }
  
  if(i == 1){
    mutations = rpois(1,lambda*ploidy)
    if(mutations > 0){
        UpdateMutationsInEachCell[[1]] = append(UpdateMutationsInEachCell[[1]],m:(m + mutations -1))
        m = m + mutations
        }

  }
  else{
    temp = 1
    for(j in WhichCellsDivided){

      mutations = rpois(1,lambda*ploidy)
      if(mutations > 0){
        UpdateMutationsInEachCell[[temp]] = append(MutationsInEachCell[[j]],m:(m+mutations-1))
        m = m + mutations
        UpdateMutationsInEachCell[[temp+1]] =  MutationsInEachCell[[j]]
      }
      else{
        UpdateMutationsInEachCell[[temp]] = MutationsInEachCell[[j]]
        UpdateMutationsInEachCell[[temp+1]] =  MutationsInEachCell[[j]]
      }
      temp = temp + 2
    }
  }
NumberOfCells = UpdateNumberOfCells
MutationsInEachCell = UpdateMutationsInEachCell

}
  
return(MutationsInEachCell)  
  
}