


Proliferation1 = function(p){
  #Proliferation1 increases a cell's probability to further divide.
  #Define new probability to divide p_d, to be given by p_d/p = 1 + s
  # where 0 <= s <= 1/p-1
  #How to choose s???
  s = rbeta(1,shape1 = 1, shape2 = 3)*(1/p-1)
  return(p*(1+s))
}


TumourGrowth = function(maxsize,beta,lambda,lambda_d,maxMuts){
#Simulate tumour until reaching "maxsize" tumour cells.
#Discrete model of tumour evolution where each cell has initially the probability of
# "beta" to divide and "1-beta" to die, so only two options. For each cell division (and therefore for each DNA replication)
# both daughter cells have the possibility to acquire new mutations. The number of passenger mutations is given by a
# Poisson distribution with rate "lambda". The number of driver mutations is also Poisson distributed with rate "lambda_d" << "lambda".
#Proliferation_function yields a recipe for how much a cell's probability of dividing increases for a random driver mutation.
#Mutations appearing at a very late stage are not expected to be detected even for driver mutations. After "maxMuts" cells, the virtual tumour grows without
#adding new mutations, since they are assumed not to be detectable anyway.

  #Begins with one tumour cell  
  NumberOfCells = 1
  
  #If driver mutations occur, cells having the same driver mutation are placed in same group:
  NumberOfCellsPerGroup = 1
  
  gen = 0 #The present generation in the branching tree.
  
  #MutationsInEachCell is what is returned. A vector of lists, where each list represent a cell, and the 
  #list contains markers (an atomic vector) for every mutation that has occured in this specific cell given by an integer.
  #Mutation 0 represents mutations already present at the first tumour cell. 
  MutationsInEachCell = 0
  
  m = 1 #Every unique passenger mutation is labeled as an integer.
  m_d = -1 #Every unique driver mutations is labeled as a negative integer.
  
  #For convenience, let the first cell be guaranteed to divide (without new driver mutations):
  NumberOfCells = 2
  NumberOfCellsPerGroup = 2
  gen = 1
  
  #Mutations in daughter cells: 
  muts1 = rpois(1,lambda)
  muts2 = rpois(1,lambda)
  #Preallocate:
  MutationsInEachCell = vector(mode = "list", length = NumberOfCells)
  if(muts1 >0){
    #add mutations:
    MutationsInEachCell[[1]] = append(0,m:(m+muts1-1))
    #update number of mutations that have occured:
    m = m + muts1
  }
  #daughter cell is just the same as mother cell:
  else{
    MutationsInEachCell[[1]] = 0
  }
  #repeat:
  if(muts2 > 0){
    MutationsInEachCell[[2]] = append(0,m:(m+muts2-1))
    m = m + muts2
  }
  else{
    MutationsInEachCell[[2]] = 0
  }
  
#From now on, let the tumour grow stochastically according to a branching process until reaching maxsize cells (or dies out).  
  
#IF TUMOUR IS NEUTRAL, OR IN OTHER WORDS lambda_d = 0:
if(lambda_d == 0){  


  while(NumberOfCells < maxsize){
 
    gen = gen + 1
  
    #Decide how many cells that divide (and therefore how many cells die)
    CellsThatDivided = rbinom(1, NumberOfCells, beta)
  
    #choose which cells divided:
    WhichCellsDivided = sample(NumberOfCells,CellsThatDivided,replace = FALSE)
  
    UpdateNumberOfCells = 2*CellsThatDivided
  
    # If tumour dies out: return(0)
    if(UpdateNumberOfCells == 0){
      return(0)
    }
  
    #preallocate vector for updating mutations in every cell (assuming this makes the code faster).
    UpdateMutationsInEachCell = vector(mode = "list", length = UpdateNumberOfCells)
  
  
    # Update mutations in each cell:
    # To avoid too much space of new mutations, stop to create new mutations after a given number of cells.
    temp = 1 #temporary variable in order to update each cell correctly.
    if(UpdateNumberOfCells < 10000000) {
      #Simulate the number of passenger mutations occuring in each of the cells
      NrOfMutationsInEachCell = rpois(UpdateNumberOfCells,lambda)
    
      for(j in WhichCellsDivided){
        #Create the two daughter cells coming from cell j:
        if(NrOfMutationsInEachCell[temp] > 0){
          UpdateMutationsInEachCell[[temp]] = append(MutationsInEachCell[[j]],m:(m+NrOfMutationsInEachCell[temp]-1))
          m = m + NrOfMutationsInEachCell[temp]
        }
        else{
          UpdateMutationsInEachCell[[temp]] = MutationsInEachCell[[j]]
        }
        if(NrOfMutationsInEachCell[temp+1] > 0){
          UpdateMutationsInEachCell[[temp+1]] = append(MutationsInEachCell[[j]],m:(m+NrOfMutationsInEachCell[temp+1]-1))
          m = m + NrOfMutationsInEachCell[temp+1]
        }
        else{
          UpdateMutationsInEachCell[[temp+1]] = MutationsInEachCell[[j]]
        }
        temp = temp + 2
      }
    }
    else{
      #UpdateMutationsInEachCell = MutationsInEachCell[rep(WhichCellsDivided,2)]
      for(j in WhichCellsDivided){
      UpdateMutationsInEachCell[[temp]] = MutationsInEachCell[[j]]
      UpdateMutationsInEachCell[[temp+1]] =  MutationsInEachCell[[j]]
      temp = temp + 2
      }
    }
    NumberOfCells = UpdateNumberOfCells
    MutationsInEachCell = UpdateMutationsInEachCell

  }
    

    return(list(MutationsInEachCell,NumberOfCells,gen))
 
  }

 ##ELSE: NON-NEUTRAL
 else{
   #Now driver mutations are allowed to occur. Whenever a cell acquires a driver mutation, create a new group that will contain all cells having this
   #spesific driver mutation (and therefore a spesific probability of dividing, p_d).
   #The first group contains all cells having p_d = beta
   MutationsInEachCellPerGroup = list(MutationsInEachCell)
   #MutationsInEachCellPerGroup[[i]] gives mutations in each cell in group i.
   #p_d[i] gives probability of dividing in group i.
   p_d = c(beta) #We assumed the first cell division gave no driver mutations for convenience.
   while(NumberOfCells < maxsize){
     gen = gen + 1

     #For each group, find how many cells divided and which cells divided:
     NrOfGroups = length(MutationsInEachCellPerGroup)
     #preallocate vectors:
     WhichCellsDividedPerGroup = vector(mode = "list",length = NrOfGroups)
     UpdateNumberOfCellsPerGroup = vector(mode = "numeric", length = NrOfGroups)
     tmp = c() #temporary variable.
     for(i in 1:NrOfGroups){
      NrOfCellsInGroupi = length(MutationsInEachCellPerGroup[[i]])
      NrOfCellsThatDividedInGroup = rbinom(1,NrOfCellsInGroupi,p_d[i])
      UpdateNumberOfCellsPerGroup[i] = 2*NrOfCellsThatDividedInGroup
      #If group i dies out:
      if(UpdateNumberOfCellsPerGroup[i] == 0){
        tmp = append(tmp,i)
      }
      #If group does not die out:
      else{
       WhichCellsDividedPerGroup[[i]] = sample(NrOfCellsInGroupi,NrOfCellsThatDividedInGroup, replace = FALSE)
      }

     }

     #Delete groups if they have died out:
     if(length(tmp)>0){
     UpdateNumberOfCellsPerGroup =  UpdateNumberOfCellsPerGroup[-tmp]
     WhichCellsDividedPerGroup = WhichCellsDividedPerGroup[-tmp]
     MutationsInEachCellPerGroup = MutationsInEachCellPerGroup[-tmp]
     p_d = p_d[-tmp]
     }
     #update total number of groups:
     NrOfGroups = length(UpdateNumberOfCellsPerGroup)
     #Update total number of cells in tumour.
     UpdateNumberOfCells = sum(UpdateNumberOfCellsPerGroup)


     # If tumour dies out: return(0)
     if(UpdateNumberOfCells == 0){
       return(0)
     }
     #After knowing which cells divided in each group, potensial mutations may now be updated:
     UpdateMutationsInEachCellPerGroup = list()
     #If total number of cells in tumour is less than maxMuts, add passenger- and driver mutations:
     if(UpdateNumberOfCells < maxMuts){

       #For each group add passenger mutations:
       for(i in 1:NrOfGroups){
         #preallocate vector for updating mutations in every cell in group i:
         UpdateMutationsInEachCellInGroupi = vector(mode = "list", length = UpdateNumberOfCellsPerGroup[i])
         #Simulate total number of mutations appearing in new generation in group i:
         #Number of passenger mutations occuring in each cell in group i:
         NrOfMutationsInEachCellInGroupi = rpois(UpdateNumberOfCellsPerGroup[i],lambda)
         temp = 1 #temporary variable.

         for(j in WhichCellsDividedPerGroup[[i]]){
           if(NrOfMutationsInEachCellInGroupi[temp] > 0){
             UpdateMutationsInEachCellInGroupi[[temp]] = append(MutationsInEachCellPerGroup[[i]][[j]],m:(m+NrOfMutationsInEachCellInGroupi[temp]-1))
             m = m + NrOfMutationsInEachCellInGroupi[temp]
           }
           else{
             UpdateMutationsInEachCellInGroupi[[temp]] = MutationsInEachCellPerGroup[[i]][[j]]
           }
           if(NrOfMutationsInEachCellInGroupi[temp+1] > 0){
             UpdateMutationsInEachCellInGroupi[[temp+1]] = append(MutationsInEachCellPerGroup[[i]][[j]],m:(m+NrOfMutationsInEachCellInGroupi[temp+1]-1))
             m = m + NrOfMutationsInEachCellInGroupi[temp+1]
           }
           else{
             UpdateMutationsInEachCellInGroupi[[temp+1]] = MutationsInEachCellPerGroup[[i]][[j]]
           }
           temp = temp + 2
         }
         UpdateMutationsInEachCellPerGroup[[i]] = UpdateMutationsInEachCellInGroupi
       }

       #Now we can check if driver mutations occur in this generation, and in that case create a new group for each cell acquiring a driver mutation.
       DriverMutations = rpois(1,lambda_d*UpdateNumberOfCells)
       if(DriverMutations > 0){
        #Distribute the number of driver mutations in each group:
        NrOfDriverMutationsInEachGroup = rmultinom(1,DriverMutations,UpdateNumberOfCellsPerGroup/UpdateNumberOfCells)
        #Find in which group(s) driver mutation(s) occur(s):
        where = which(NrOfDriverMutationsInEachGroup>0)
        #Decide which cell(s) that acquire(s) driver mutation(s) in each group:
        for(i in where){
          #Pick cells randomly assuming a new-born cell cannot acquire more than one driver mutation for convenience (therefore replace = FALSE):
          WhichCellsAcquireDriverMutationInGroupi = sample(UpdateNumberOfCellsPerGroup[i],NrOfDriverMutationsInEachGroup[i],replace = FALSE)
          #Update NumberOfCellsPerGroup:
          UpdateNumberOfCellsPerGroup[i] = UpdateNumberOfCellsPerGroup[i] - length(WhichCellsAcquireDriverMutationInGroupi)
          for(j in WhichCellsAcquireDriverMutationInGroupi){
            #Create new group:
            UpdateMutationsInEachCellPerGroup[[length(UpdateMutationsInEachCellPerGroup)+1]] = list(append(UpdateMutationsInEachCellPerGroup[[i]][[j]],m_d:(m_d-NrOfDriverMutationsInEachGroup[i]+1)))
            m_d = m_d - NrOfDriverMutationsInEachGroup[i]
            #Add any proliferation factor
            p_d[length(p_d)+1] = Proliferation1(p_d[i])
            #Update NumberOfCellsPerGroup:
            UpdateNumberOfCellsPerGroup[length(UpdateNumberOfCellsPerGroup)+1] = 1
          }
          #delete cells from former group:
          UpdateMutationsInEachCellPerGroup[[i]] = UpdateMutationsInEachCellPerGroup[[i]][-WhichCellsAcquireDriverMutationInGroupi]
        }

       }

     }
     #Else, do not add any mutations, as they are likely not to be traceable anyway.
     else{
       for(i in 1:NrOfGroups){
         #preallocate vector for updating mutations in every cell in group i:
         UpdateMutationsInEachCellInGroupi = vector(mode = "list", length = UpdateNumberOfCellsPerGroup[i])
         temp = 1
         for(j in WhichCellsDividedPerGroup[[i]]){
          UpdateMutationsInEachCellInGroupi[[temp]] = MutationsInEachCellPerGroup[[i]][[j]]
          UpdateMutationsInEachCellInGroupi[[temp+1]] = MutationsInEachCellPerGroup[[i]][[j]]
          temp = temp + 2
         }
         UpdateMutationsInEachCellPerGroup[[i]] = UpdateMutationsInEachCellInGroupi
       }
     }
     #Update mutations:
     MutationsInEachCellPerGroup = UpdateMutationsInEachCellPerGroup
     #Update total number of cells
     NumberOfCells = UpdateNumberOfCells
     #NumberOfCellsPerGroup = UpdateNumberOfCellsPerGroup
   }
   return(list(MutationsInEachCellPerGroup,NumberOfCells,NumberOfCellsPerGroup,gen,p_d))
 }
  
}

TakeBiopsy = function(TumourSample, S, ReadDepth, ploidy = 2){
  #Function has input from a tumour, and take a biopsy consisting of S cells.
  #The tumour cells are assumed to be well-mixed in the tumour.
  #The average ploidy is assumed to be near 2.
  
  #TumourSample may be a list of vectors of lists (for non-neutral tumours) or a vector of lists.
  TotalNumberOfCells = TumourSample[[2]]
  #Take a biopsy of S cells:
  #If tumour is neutral:
  if(length(TumourSample)==3){
    #Sample S cells from tumour. Assuming the cells are well-mixed all cells are equally likely to be chosen: 
    TheCellsChosen = sample(TotalNumberOfCells,S,replace = FALSE)
    TheSample = TumourSample[[1]][TheCellsChosen] 
  }
  #else, tumour is non-neutral:
  else{
    #Count how many cells are sampled from each group according to a multinomial distribution:
    NrOfCellsChosenInEachGroup = rmultinom(1,S, TumourSample[[3]]/TumourSample[[2]])
    #Which groups have one or more cells represented:
    where = which(NrOfCellsChosenInEachGroup>0)
    #Preallocate:
    TheSample = vector(mode = "list", length = S)
    temp = 1
    for(i in where){
      CellsChosenInGroupi = sample(TumourSample[[3]][i],NrOfCellsChosenInEachGroup[i],replace = FALSE)
      TheSample[temp:(temp + NrOfCellsChosenInEachGroup[i]-1)] = TumourSample[[1]][[i]][CellsChosenInGroupi]
      temp = temp + NrOfCellsChosenInEachGroup[i]
    }
  
  }
  
  #Get all mutation IDs in one single atomic vector:
  AllMutations = unlist(TheSample)
  #Get VAF of all somatic mutations using R-function table:
  RealVAF = table(AllMutations)/(S*ploidy)
  #Frequencies that are too low are not detectable and have too much uncertainty:
  RealVAF = RealVAF[RealVAF >= 0.01]
  
  #Now, the OBSERVED frequencies of each somatic mutation is estimated assuming a constant read depth:
  ObservedVAF = vector(mode = "numeric",length = length(RealVAF))
  for(i in 1:length(RealVAF)){
    NrOfMutationsRecored = rbinom(1,ReadDepth,as.numeric(RealVAF[i]))
    ObservedVAF[i] =NrOfMutationsRecored/ReadDepth
  }
  
  
  #hist(ObservedVAF, freq = 0, breaks = 50)
  return(ObservedVAF)
}


library(foreach)
library(doParallel)
cl = makeCluster(10)
registerDoParallel(cl)
foreach(i = 1:30) %dopar% {
  
  tumourD = TumourGrowth(100000000,0.55,1.5,10^-8,10000000)
  
  if(object.size(tumourN) > 70){
    biopsy = TakeBiopsy(tumourD,1000000,100)
    nameN = paste("/home/shomeb/p/paalvj/Documents/Masteroppgave/","tumourD",i,".RData",sep = "")
    save(biopsy, file = nameN)
  }
  
}
stopCluster(cl)




