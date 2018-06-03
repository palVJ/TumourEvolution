Proliferation = function(p,s){
  #Proliferation increases a cell's probability to further divide.
  #Define new probability to divide p_d, to be given by p_d/p = 1 + s
  # where 0 <= s <= 1/p-1
  
  p_d = p*(1+s)
  if(p_d <= 1){
    return(p_d)
  }
  #Probability of cell division can't be higher than 1. 
  else{
    return(p)
  }
}


TumourGrowth = function(maxgen,beta,lambda,lambda_d,maxMuts,s){
#Simulate tumour until reaching generation maxgen.
#Discrete branching tree model of tumour evolution where each cell has initially the probability of
# "beta" to divide and "1-beta" to die, so only two options. For each cell division
# both daughter cells have the possibility to acquire new mutations. The number of passenger mutations is given by a
# Poisson distribution with rate "lambda". The number of driver mutations is also Poisson distributed with rate "lambda_d" << "lambda".
#Proliferation function gives a recipe for how much a cell's probability of dividing increases for a random driver mutation.
#s gives the proliferation factor by which the cell division probability increases by p(1+s). Se proliferation function above. 
#Mutations appearing at a very late stage are not expected to be detected even for driver mutations. After "maxMuts" cells, the virtual tumour grows without
#adding new mutations, since the mutations are assumed not to be traceable anyway.
  
#The function returns the mutations for all cells, the last generation in the branching tree, the total number of cells, the number of cells in each 
#group having the same probability of dividing, and lastly the probability of dividing for each group. 
#If lamdda_d = 0, this means the tumour grows neutrally. 

  #Begins with one tumour cell  
  NumberOfCells = 1
  
  #If driver mutations occur, cells having the same driver mutation are placed in same group. Now there is only one cell in one group:
  NumberOfCellsPerGroup = 1
  
  gen = 0 #The present generation in the branching tree. The initial generation is by definition 0 (consisting of the first tumour cell)
  
  #MutationsInEachCell is what is returned. For neutral evolution: A vector of lists, where each list represent a cell, and the 
  #list contains an atomic vector consisting of every mutation that has occured in this specific cell given by an integer.
  #For non-neutral evolution: A list of vectors, where each vector consists of lists. Each vector represent a group of cells having the same fitness.
  #Each list in a vector represents a cell in this group. The list consists of all mutations this cell in this group has. 
  #Mutation 0 represents mutations already present at the first tumour cell. 
  MutationsInEachCell = 0
  m = 1 #Every unique passenger mutation is labeled as an integer.
  m_d = -1 #Every unique driver mutation is labeled as a negative integer.
  
  #For convenience, let the first cell be guaranteed to divide without new driver mutations:
  NumberOfCells = 2
  NumberOfCellsPerGroup = 2
  gen = 1
  
  #Mutations in daughter cells: 
  muts1 = rpois(1,lambda)
  muts2 = rpois(1,lambda)
  #Preallocate:
  MutationsInEachCell = vector(mode = "list", length = NumberOfCells)
  InWhichGenerationMutationsOccurred = vector(mode = "list", length = NumberOfCells)
  if(muts1 >0){
    #add mutations:
    MutationsInEachCell[[1]] = c(0,m:(m+muts1-1))
    InWhichGenerationMutationsOccurred[[1]] = c(0,rep(gen,muts1))
    #update number of mutations that have occured:
    m = m + muts1
  }
  #daughter cell is just the same as mother cell:
   else{
     MutationsInEachCell[[1]] = 0
     InWhichGenerationMutationsOccurred[[1]] = 0
   }
  #repeat:
  if(muts2 > 0){
    MutationsInEachCell[[2]] = c(0,m:(m+muts2-1))
    InWhichGenerationMutationsOccurred[[2]] = c(0,rep(gen,muts2))
    m = m + muts2
  }
   else{
     MutationsInEachCell[[2]] = 0
     InWhichGenerationMutationsOccurred[[2]] = 0
   }
  
#From now on, let the tumour grow stochastically according to a branching process until reaching maxsize cells (or dies out).  
  
#IF TUMOUR IS NEUTRAL, OR IN OTHER WORDS lambda_d = 0:

#if(lambda_d == 0){

  while(gen <= maxgen){
 
    
  
    #Decide how many cells that divide in present generation (the other cells will eventually die):
    CellsThatDivided = rbinom(1, NumberOfCells, beta)
  
    #choose which cells divided:
    WhichCellsDivided = sample(NumberOfCells,CellsThatDivided,replace = FALSE)
  
    UpdateNumberOfCells = 2*CellsThatDivided
  
    # If tumour dies out: return(0)
    if(UpdateNumberOfCells == 0){
      return(0)
    }
    
    #Since tumour has not died out, the generation has increased by one
    gen = gen + 1
    
    #preallocate vector for updating mutations in every cell (assuming this makes the code faster).
    UpdateMutationsInEachCell = vector(mode = "list", length = UpdateNumberOfCells)
    UpdateInWhichGenerationMutationsOccurred = vector(mode = "list", length = UpdateNumberOfCells)
  
    # Update mutations in each cell:
    # To avoid too much space of new mutations, stop to create new mutations after a given generation given by maxMuts
    temp = 1 #temporary variable in order to update each cell correctly.
    if(gen <= maxMuts) {
      #Simulate the number of passenger mutations occuring in each of the cells
      NrOfMutationsInEachCell = rpois(UpdateNumberOfCells,lambda)
      
      #Update mutations in each cell:
      for(j in WhichCellsDivided){
        #Create the two daughter cells coming from cell j:
        if(NrOfMutationsInEachCell[temp] > 0){
          UpdateMutationsInEachCell[[temp]] = c(MutationsInEachCell[[j]],m:(m+NrOfMutationsInEachCell[temp]-1))
          UpdateInWhichGenerationMutationsOccurred[[temp]] = c(InWhichGenerationMutationsOccurred[[j]],rep(gen,NrOfMutationsInEachCell[temp]))
          m = m + NrOfMutationsInEachCell[temp]
        }
        else{
          UpdateMutationsInEachCell[[temp]] = MutationsInEachCell[[j]]
          UpdateInWhichGenerationMutationsOccurred[[temp]] = InWhichGenerationMutationsOccurred[[j]]
        }
        if(NrOfMutationsInEachCell[temp+1] > 0){
          UpdateMutationsInEachCell[[temp+1]] = c(MutationsInEachCell[[j]],m:(m+NrOfMutationsInEachCell[temp+1]-1))
          UpdateInWhichGenerationMutationsOccurred[[temp+1]] = c(InWhichGenerationMutationsOccurred[[j]],rep(gen,NrOfMutationsInEachCell[temp+1]))
          m = m + NrOfMutationsInEachCell[temp+1]
        }
        else{
          UpdateMutationsInEachCell[[temp+1]] = MutationsInEachCell[[j]]
          UpdateInWhichGenerationMutationsOccurred[[temp+1]] = InWhichGenerationMutationsOccurred[[j]]
        }
        temp = temp + 2
      }
     }
    #else, there is no need to add mutations as they will not be detectable. 
    else{
      #All daughter cells inherit exactly the same mutations as their parents:
      UpdateMutationsInEachCell = MutationsInEachCell[rep(WhichCellsDivided,2)]
      UpdateInWhichGenerationMutationsOccurred = InWhichGenerationMutationsOccurred[rep(WhichCellsDivided,2)]
      # for(j in WhichCellsDivided){
      # UpdateMutationsInEachCell[[temp]] = MutationsInEachCell[[j]]
      # UpdateMutationsInEachCell[[temp+1]] =  MutationsInEachCell[[j]]
      # temp = temp + 2
      # }
    }
    #Update parameters
    NumberOfCells = UpdateNumberOfCells
    MutationsInEachCell = UpdateMutationsInEachCell
    InWhichGenerationMutationsOccurred = UpdateInWhichGenerationMutationsOccurred
    }

    return(list(MutationsInEachCell,InWhichGenerationMutationsOccurred,NumberOfCells,gen))
 
  }

  
 # ##ELSE: NON-NEUTRAL
 # else{
 #   #Now driver mutations are allowed to occur. Whenever a cell acquires a driver mutation, create a new group that will contain all cells having this
 #   #spesific driver mutation (and therefore a spesific probability of dividing, p_d).
 #   #The first group contains all cells having p_d = beta
 #   
 #   #Create a list containing the initial vector MutationsInEachCell consisting of the two cells.
 #   MutationsInEachCellPerGroup = list(MutationsInEachCell)
 #   #MutationsInEachCellPerGroup[[i]] gives mutations in each cell in group i.
 #   #p_d[i] gives probability of dividing in group i.
 #   p_d = c(beta) #We assumed the first cell division gave no driver mutations for convenience.
 #   while(NumberOfCells < maxsize){
 #     
 # 
 #     #Find how many groups there are:
 #     NrOfGroups = length(MutationsInEachCellPerGroup)
 #     
 #     #For each group, find how many cells divided and which cells divided:
 #     #preallocate vectors:
 #     WhichCellsDividedPerGroup = vector(mode = "list",length = NrOfGroups)
 #     UpdateNumberOfCellsPerGroup = vector(mode = "numeric", length = NrOfGroups)
 #     tmp = c() #temporary variable.
 #     for(i in 1:NrOfGroups){
 #      NrOfCellsInGroupi = length(MutationsInEachCellPerGroup[[i]])
 #      NrOfCellsThatDividedInGroup = rbinom(1,NrOfCellsInGroupi,p_d[i])
 #      UpdateNumberOfCellsPerGroup[i] = 2*NrOfCellsThatDividedInGroup
 #      #If group i dies out:
 #      if(UpdateNumberOfCellsPerGroup[i] == 0){
 #        tmp = c(tmp,i)
 #      }
 #      #If group i does not die out:
 #      else{
 #       WhichCellsDividedPerGroup[[i]] = sample(NrOfCellsInGroupi,NrOfCellsThatDividedInGroup, replace = FALSE)
 #      }
 # 
 #     }
 # 
 #     #Delete groups if they have died out:
 #     if(length(tmp)>0){
 #     UpdateNumberOfCellsPerGroup =  UpdateNumberOfCellsPerGroup[-tmp]
 #     WhichCellsDividedPerGroup = WhichCellsDividedPerGroup[-tmp]
 #     MutationsInEachCellPerGroup = MutationsInEachCellPerGroup[-tmp]
 #     p_d = p_d[-tmp]
 #     }
 #     #update total number of groups:
 #     NrOfGroups = length(UpdateNumberOfCellsPerGroup)
 #     #Update total number of cells in tumour.
 #     UpdateNumberOfCells = sum(UpdateNumberOfCellsPerGroup)
 # 
 # 
 #     # If tumour dies out: return(0)
 #     if(UpdateNumberOfCells == 0){
 #       return(0)
 #     }
 #     
 #     #Since tumour did not die out, increase generation by one:
 #     gen = gen + 1
 #     
 #     #After knowing which cells divided in each group, potensial new mutations may now be added:
 #     UpdateMutationsInEachCellPerGroup = list()
 #     #If total number of cells in tumour is less than maxMuts, add passenger- and driver mutations:
 #     if(UpdateNumberOfCells < maxMuts){
 # 
 #       #For each group add passenger mutations:
 #       for(i in 1:NrOfGroups){
 #         #preallocate vector for updating mutations in every cell in group i:
 #         UpdateMutationsInEachCellInGroupi = vector(mode = "list", length = UpdateNumberOfCellsPerGroup[i])
 #         #Simulate total number of mutations appearing in new generation in group i:
 #         #Number of passenger mutations occuring in each cell in group i:
 #         NrOfMutationsInEachCellInGroupi = rpois(UpdateNumberOfCellsPerGroup[i],lambda)
 #         temp = 1 #temporary variable.
 #          
 #         #For each cell j in group i that divided in the former generation:
 #         for(j in WhichCellsDividedPerGroup[[i]]){
 #           #if daughter cell number 1. acquires passenger mutations:
 #           if(NrOfMutationsInEachCellInGroupi[temp] > 0){
 #             #create daughter cell having the same mutations as cell j in addition to new mutations.
 #             UpdateMutationsInEachCellInGroupi[[temp]] = c(MutationsInEachCellPerGroup[[i]][[j]],m:(m+NrOfMutationsInEachCellInGroupi[temp]-1))
 #             m = m + NrOfMutationsInEachCellInGroupi[temp]
 #           }
 #           #Else, create daughter cell having exactly the same mutations as parent. 
 #           else{
 #             UpdateMutationsInEachCellInGroupi[[temp]] = MutationsInEachCellPerGroup[[i]][[j]]
 #           }
 #           #Do exactly the same for daughter cell number 2: 
 #           if(NrOfMutationsInEachCellInGroupi[temp+1] > 0){
 #             UpdateMutationsInEachCellInGroupi[[temp+1]] = c(MutationsInEachCellPerGroup[[i]][[j]],m:(m+NrOfMutationsInEachCellInGroupi[temp+1]-1))
 #             m = m + NrOfMutationsInEachCellInGroupi[temp+1]
 #           }
 #           else{
 #             UpdateMutationsInEachCellInGroupi[[temp+1]] = MutationsInEachCellPerGroup[[i]][[j]]
 #           }
 #           temp = temp + 2
 #         }
 #         UpdateMutationsInEachCellPerGroup[[i]] = UpdateMutationsInEachCellInGroupi
 #       }
 # 
 #       #Now we can check if driver mutations occur in this new generation, and in that case create a new group for each cell acquiring a driver mutation.
 #       #Number of driver mutations occuring in each group:
 #       NrOfDriverMutationsInEachGroup = rpois(NrOfGroups,lambda_d*UpdateNumberOfCellsPerGroup)
 #       #If driver mutations occurred:
 #       if(sum(NrOfDriverMutationsInEachGroup) > 0){
 #        #check in which groups driver mutations occurred:
 #        where = which(NrOfDriverMutationsInEachGroup>0)
 #        #Decide which cell(s) that acquire(s) driver mutation(s) for each group i:
 #        for(i in where){
 #          #Pick cells randomly to acquire driver mutations:
 #          #Find out which cells that acquire driver mutations:
 #          WhichCellsAcquireDriverMutationInGroupi = unique(sample(UpdateNumberOfCellsPerGroup[i],NrOfDriverMutationsInEachGroup[i],replace = TRUE))
 #          #Comment on code above: replace = TRUE in case a cell gets more than one driver mutation (very very unlikely, but most be accounted for, or else code can fail)
 #          #unique, because for convenience, we don't differ between cells having one driver mutation or more than one driver mutation. 
 #          
 #          #Update Number of cells in group i (move cells that aquired driver mutations to seperate groups):
 #          UpdateNumberOfCellsPerGroup[i] = UpdateNumberOfCellsPerGroup[i] - length(WhichCellsAcquireDriverMutationInGroupi)
 #          
 #          #Add a new group for each driver mutation:
 #          for(j in 1:length(WhichCellsAcquireDriverMutationInGroupi)){
 #            UpdateMutationsInEachCellPerGroup[[length(UpdateMutationsInEachCellPerGroup)+1]] = 
 #            list(c(UpdateMutationsInEachCellPerGroup[[i]][[WhichCellsAcquireDriverMutationInGroupi[j]]],m_d))
 #            m_d = m_d -1
 #            #Add any proliferation factor increasing the cell's probability to divide:
 #            p_d[length(p_d)+1] = Proliferation(p_d[i],s)
 #            #Update NumberOfCellsPerGroup for this group now consisting of only one cell:
 #            UpdateNumberOfCellsPerGroup[length(UpdateNumberOfCellsPerGroup)+1] = 1
 #          }
 #          
 #          #delete cells from former group as they are now removed to a new group:
 #          UpdateMutationsInEachCellPerGroup[[i]] = UpdateMutationsInEachCellPerGroup[[i]][-WhichCellsAcquireDriverMutationInGroupi]
 #          
 #        }
 # 
 #       }
 # 
 #     }
 #     
 #     #Else, if NumberOfCells > MaxMuts, do not add any mutations, as they are likely not to be traceable anyway.
 #     else{
 #       for(i in 1:NrOfGroups){
 #         
 #         UpdateMutationsInEachCellInGroupi = MutationsInEachCellPerGroup[[i]][rep(WhichCellsDividedPerGroup[[i]],2)]
 #         # #preallocate vector for updating mutations in every cell in group i:
 #         # UpdateMutationsInEachCellInGroupi = vector(mode = "list", length = UpdateNumberOfCellsPerGroup[i])
 #         # temp = 1
 #         # for(j in WhichCellsDividedPerGroup[[i]]){
 #         #  #Both daughter cells inherit exactly the same mutations as their parent:
 #         #  UpdateMutationsInEachCellInGroupi[[temp]] = MutationsInEachCellPerGroup[[i]][[j]]
 #         #  UpdateMutationsInEachCellInGroupi[[temp+1]] = MutationsInEachCellPerGroup[[i]][[j]]
 #         #  temp = temp + 2
 #         # }
 #         UpdateMutationsInEachCellPerGroup[[i]] = UpdateMutationsInEachCellInGroupi
 #       }
 #     }
 #     #Update mutations:
 #     MutationsInEachCellPerGroup = UpdateMutationsInEachCellPerGroup
 #     #Update total number of cells
 #     NumberOfCells = UpdateNumberOfCells
 #     NumberOfCellsPerGroup = UpdateNumberOfCellsPerGroup
 #   }
 #   
 #   return(list(MutationsInEachCellPerGroup,NumberOfCells,NumberOfCellsPerGroup,gen,p_d))
 # }
 #}
  
  RealDistributionOfVAF = function(TumourSample,gens,ploidy = 2){
    
    TotalNumberOfCells = TumourSample[[3]]
    TheSample = TumourSample[[1]]
    CorGen = TumourSample[[2]]
    #Get all mutation IDs in one single atomic vector:
    AllMutations = unlist(TheSample)
    #Get VAF of all somatic mutations using R-function table:
    RealVAF = table(AllMutations)/(ploidy*TotalNumberOfCells)
    #Get the corresponding generation where each mutation occurred:
    GenerationForEachVaf = unlist(CorGen)
    Vafs = vector(mode = "list", length = gens)
    for(i in 1:gens){
      #Find which mutations belong to the same generation (some mutations are equal as they belong to the same parent possesing this mutation)
      where = which(GenerationForEachVaf == i)
      #Find the unique mutations appearing in generation i
      muts = as.character(unique(AllMutations[where]))
      #Store the corresponding VAFs of the somatic mutations appearing in generation i.
      Vafs[[i]] = as.vector(RealVAF[muts])
    }
    return(Vafs)
  }
  
  #COMPUTE OBSERVED M(k) from tumor:
  RealMk = function(TumourSample,gens){
    #tumor contains cells and mutations in each cell in addition to which generation each surviving mutation appeared.
    #gens is the number of generations to look at. 
    TotalNumberOfCells = TumourSample[[3]]
    TheSample = TumourSample[[1]]
    CorGen = TumourSample[[2]]
    #Get all mutation IDs in one single atomic vector:
    AllMutations = unlist(TheSample)
    #Get the corresponding generation where each mutation occurred:
    GenerationForEachVaf = unlist(CorGen)
    NrOfMutsAppearedInGenk = vector(mode = "numeric", length = gens)
    for(k in 1:gens){
      #Find which mutations belong to the same generation (some mutations are equal as they belong to the same parent possesing this mutation)
      where = which(GenerationForEachVaf == k)
      #Find the unique number of mutations appearing in generation i
      NrOfMutsAppearedInGenk[k] = length(unique(AllMutations[where]))
    }
    Mk = cumsum(NrOfMutsAppearedInGenk)
    return(Mk)
  }
  

TakeResection = function(TumourSample, S, par_mean,par_size,purity, ploidy = 2){
  #Function has input "TumourSample" from a virtual tumour created from TumourGrowth(). Take a biopsy consisting of S cells.
  #The tumour cells are assumed to be well-mixed in the tumour (a strong assumption!!).
  #purity is the amount of cells in S that is tumour cells
  #par_mean and size is the parameters of the negative binomial distribution (mean and size). 
  
  
  #TumourSample may be a list of vectors of lists (for non-neutral tumours) or a vector of lists (neutral).
  #Find total number of cells in virtual tumour:
  TotalNumberOfCells = TumourSample[[2]]
  #Take a biopsy of S cells:
  #If tumour is neutral:
  if(length(TumourSample)==3){
    #Sample S cells from tumour of which S*purity are tumour cells. Assuming the cells are well-mixed all cells are equally likely to be chosen: 
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
  RealVAF = (table(AllMutations)/(ploidy*S))*purity
  #Frequencies that are too low have very little chance of beeing seen, especially in the region of interest [0.01,0.25].
  #Therefore, for simplicity discard too low frequencies as they will not influence the great picture anyway:
  RealVAF = RealVAF[RealVAF >= 0.001]
  
  #Now, the OBSERVED frequencies of each somatic mutation is estimated assuming a poisson distributed read depth:
  ObservedVAF = vector(mode = "numeric",length = length(RealVAF))
  #Simulate read depth in each position:
  ReadDepths = rnbinom(n = length(RealVAF), size = par_size, mu = par_mean )
  for(i in 1:length(RealVAF)){
   
    if(ReadDepths[i] >= 10){
    NrOfMutationsRecorded = rbinom(1,ReadDepths[i],as.numeric(RealVAF[i]))
    
      if(NrOfMutationsRecorded >= 3){
        ObservedVAF[i] = NrOfMutationsRecorded/ReadDepths[i]
    }
    
    }
    
  }
  
  #Look only at the frequency region at interest:
  #Lower bound:
  ObservedVAF = ObservedVAF[ObservedVAF >= 0.01]
  
  
  
  return(list(RealVAF,ObservedVAF))
}

rplbounded = function(n,alpha, xmin, xmax){
  #Simulates n times from a bounded power law distribution
  #with exponent rate alpha, and interval [xmin,xmax].
  
  sims = vector(mode = "numeric", length = n)
  unis = runif(n)
  for(i in 1:n){
    sims[i] = (unis[i]*(xmax^(1-alpha)-xmin^(1-alpha))+xmin^(1-alpha))^(1/(1-alpha))
  }
  return(sims)
}

VirtualSequenceAnalysis = function(n,alpha,xmin,xmax,par_mean,par_size,purity){
  
  #Simulate the number of somatic mutations given in the tissue sample in addition to the real VAF
  #for each somatic mutation in the region [xmin,xmax]:
  #VAFmutations = rplbounded(n,alpha,xmin,xmax)
  VAFmutations = rexp(n, rate = 1)
  VAF = purity*VAFmutations
  ObservedVAF = vector(mode = "numeric", length = n)
  Coverages = rnbinom(n,size = par_size, mu = par_mean)
  for(i in 1:n){
    if(Coverages[i] >= 8){
      observed = rbinom(1,Coverages[i],VAF[i])
      if(observed >= 3){
        ObservedVAF[i] = (observed/Coverages[i])*purity
      }
    }
  }
  ObservedVAF = ObservedVAF[ObservedVAF >= 0.01]
  return(ObservedVAF)
}

DistributionOfPublicMutations = function(ID,readsFile, purity, sims,l){
  
  #Get read depth data.
  name = paste(readsFile,ID, sep = "/")
  name = paste(name, ".RData",sep = "")
  readsData = load(file = name)
  reads = get(readsData)
  #Compute MLE-estimators:
  library(MASS)
  par_estimators = fitdistr(x = reads, densfun = "negative binomial")
  size_par = as.numeric(par_estimators$estimate[1])
  mu_par = as.numeric(par_estimators$estimate[2])
    
  #Generate sims reads:
  r = rnbinom(n = sims, size = size_par, mu = mu_par)
  #For each read, generate an observed VAF-value:
  VAF_observed = rbinom(n = sims, size = r, prob = 0.5*purity)/r
  p = length(which(VAF_observed < l))/sims
  hist(VAF_observed, breaks = 100, xlim = c(0,0.6))
  return(p)
}

DistributionOfExtinction = function(beta, maxgen){
#Count number of cells in each generation until extinction:
#population begins with one cell:
pop = 1
gen = 1
NrOfCellsInEachGen = vector(mode = "numeric", length = maxgen)
while(pop>0 && gen <= maxgen){
#Decide how many cells that divide in present generation (the other cells will eventually die):
CellsThatDivided = rbinom(1,pop, beta)

pop = 2*CellsThatDivided
# If tumour dies out: return(0)
if(pop == 0){
  return(NrOfCellsInEachGen)
}
else{
  NrOfCellsInEachGen[gen] = pop
}
gen = gen + 1
}
  #extinction not guarantied, return 0:
  return(0)
}

DistributionOfSurivingPopulations = function(beta,maxgen){
  #Count number of cells in each generation until maxgen
  #population begins with one cell:
  pop = 1
  gen = 1
  NrOfCellsInEachGen = vector(mode = "numeric", length = maxgen)
  while(gen <= maxgen){
    #Decide how many cells that divide in present generation (the other cells will eventually die):
    CellsThatDivided = rbinom(1,pop, beta)
  
    pop = 2*CellsThatDivided
    # If tumour dies out: return(0)
    if(pop == 0){
      return(0)
    }
    else{
      NrOfCellsInEachGen[gen] = pop
    }
    gen = gen + 1
  }
  return(NrOfCellsInEachGen)
}

OrderNumericInRightPlace = function(x,n) {
  # x is sorted vector. n is number to be placed in right position.
  a = x
  while (length(a) > 1) {
    
    if (n >= a[length(a) / 2]) {
      a = a[1:(length(a)/2)]
    }
    else {
      a = a[((length(a)/2)+1):length(a)]
    }
    
  }

  if (n > a) {
    position = which(x == a)-1
  }
  else {
    position = which(x == a)
  }
  # Returns between which elements in vector x one must place n.
  return(position)
}





#Decide parameters here:
beta = 0.55
lambda = 1.5
lambda_d = 0
s = 0
maxgen = 160
maxMuts = 50
ploidy = 2
par_mean = c(80,90,100,110,120,130,190)
par_size = c(2.5,2.5,2.5,2.5,2.5,2.5,2.5)
purity = c(0.7,0.7,0.7,0.7,0.7,0.7,0.7)


###################################################################
#CREATE A TUMOR AND COMPUTE M(E[F_k])
  k = 1:40
  Efk = 1/(2*(2*beta)^k)
  x = 1/Efk
  mEfk = (lambda/(2-1/beta))*(1/Efk-2)
  t = TumourGrowth(maxgen,beta,lambda,lambda_d,maxMuts, s)
    while(object.size(t)<70){
      t= TumourGrowth(maxgen,beta,lambda,lambda_d,maxMuts, s )
    }
  TotalNumberOfCells = t[[3]]
  TheSample = t[[1]]
  #Get all mutation IDs in one single atomic vector:
  AllMutations = unlist(TheSample)
  #Get VAF of all somatic mutations using R-function table:
  RealVAF = table(AllMutations)/(ploidy*TotalNumberOfCells)
  RealVAF = RealVAF[RealVAF >= 0.01]
  tab = rev(table(RealVAF))
  Mf = as.vector(cumsum(tab))
  recordedvafs = as.numeric(names(tab))
  MEfk = vector(mode = "numeric", length = length(Efk))
  for(i in 1:length(MEfk)){
    position = OrderNumericInRightPlace(recordedvafs,Efk[i])
    MEfk[i] = Mf[position]
  }
  Mflist = list(x,mEfk, MEfk)
  save(Mflist, file = paste("/home/shomeb/p/paalvj/Documents/Masteroppgave/","Mflist2",beta,lambda,".RData",sep = "_"))
###################################################################

#COMPARE M(f) in real data vs. observed data with given read depth distribution:
# 
#  library(foreach)
#  library(doParallel)
#  #Decide number of clusters here:
#  cl = makeCluster(3)
#  registerDoParallel(cl)
#  
#  foreach(i = 1:7) %dopar% {
# o = TakeResection(t,S,par_mean[i],par_size[i],purity[i])
# save(o, file = paste("/home/shomeb/p/paalvj/Documents/Masteroppgave/","observedVAFVsRealVAFOnlyChangemean2",par_mean[i],par_size[i],purity[i],".RData",sep = "_"))
#  }
#  stopCluster(cl)
########################################################################################
  

#ESTIMATE distribution of VAF here when neutral evolution (only looking at subclonal mutations):
# NrOfTumours = 200
#  library(foreach)
#  library(doParallel)
#  #Decide number of clusters here:
#  cl = makeCluster(5)
#  registerDoParallel(cl)
# foreach(i = 1:NrOfTumours) %dopar% {
#   t = TumourGrowth(maxgen,beta,lambda,lambda_d,maxMuts, s)
#   while(object.size(t)<70){
#     t= TumourGrowth(maxgen,beta,lambda,lambda_d,maxMuts, s )
#   }
#   r = RealDistributionOfVAF(TumourSample = t,gens = maxMuts)
#   save(r, file =  paste("/home/shomeb/p/paalvj/Documents/Masteroppgave/EstVAFDist5/",beta,lambda,i,".RData",sep= "_"))
# }
##############################################################################################
  
# #ESTIMATE EXPECTED GROWTH FOR SURVIVING TUMORS
# NrOfTumours = 2000
# library(foreach)
# library(doParallel)
# #Decide number of clusters here:
# cl = makeCluster(5)
# registerDoParallel(cl)
# foreach(i = 1:NrOfTumours) %dopar% {
#   sur = DistributionOfSurivingPopulations(beta,maxgen)
#   if(object.size(sur) > 100){
#     save(sur, file =  paste("/home/shomeb/p/paalvj/Documents/Masteroppgave/SurvivingDistribution/",beta,lambda,i,".RData",sep= "_"))
#   }
# }

#########################################################################################
# #FIND M(k) based on 200 tumors
# NrOfTumours = 100
# gens = 40
#  library(foreach)
#  library(doParallel)
#  #Decide number of clusters here:
#  cl = makeCluster(5)
#  registerDoParallel(cl)
# foreach(i = 1:NrOfTumours) %dopar% {
#   t = TumourGrowth(maxgen,beta,lambda,lambda_d,maxMuts, s)
#   while(object.size(t)<70){
#     t= TumourGrowth(maxgen,beta,lambda,lambda_d,maxMuts, s )
#   }
#   Mk = RealMk(t,gens)
#   save(Mk, file =  paste("/home/shomeb/p/paalvj/Documents/Masteroppgave/Mk3/",beta,lambda,300+i,".RData",sep = "_"))
# }
#######################################################################################
  
#ESTIMATE expected number of cells in each generation conditioned on extinction:
# maxgen = 120
# extinctMatrix = matrix(nrow = 5000000, ncol = maxgen)
# WhenDied = vector(mode = "numeric", length = 120)
# for(i in 1:5000000){
#   d = DistributionOfExtinction(beta,maxgen)
#   if(length(d)>1){
#     d = DistributionOfExtinction(beta,maxgen)
#     when = min(which(d == 0))
#     WhenDied[when] =  WhenDied[when] + 1
#     extinctMatrix[i,] = d
#   }
# 
# }
# ex = colMeans(extinctMatrix, na.rm = TRUE)
# save(ex,file = "/home/shomeb/p/paalvj/Documents/Masteroppgave/ex.RData")
# save(WhenDied, file = "/home/shomeb/p/paalvj/Documents/Masteroppgave/WhenDied.RData")

#######################################################################################

#CHECH DISTRIBUTION OF PUBLIC MUTATIONS:
# pm = DistributionOfPublicMutations(ID = "TCGA-NH-A5IV-01A-42D-A36X-10", readsFile = "/home/shomeb/p/paalvj/Documents/Masteroppgave/Read_Depths",
#                                   purity = 0.8,sims = 1000000, l = 0.1)



#v = VirtualSequenceAnalysis(200,2,0.01,0.5,200,2.4,0.72)
# tumour = TumourGrowth(maxsize,beta,lambda,lambda_d,maxMuts,s)
# 
#   while(object.size(tumour)<70){
#     tumour = TumourGrowth(maxsize,beta,lambda,lambda_d,maxMuts, s )
#   }
# biopsy = TakeBiopsy(tumour,10000000,100,2,1)
# save(biopsy, file = paste("/home/shomeb/p/paalvj/Documents/Masteroppgave/","driverbiopsy2","100mTo1m",beta,lambda,lambda_d,s,".RData",sep = "_"))
# # #Simulate bunches of tumours in parallell:
# library(foreach)
# library(doParallel)
# #Decide number of clusters here:
# cl = makeCluster(3)
# registerDoParallel(cl)
# 
# foreach(i = 1:11) %dopar% {
# 
# #Call TumourGrowth() until a tumour has grown to desired size
#   tumour = TumourGrowth(maxsize,beta,lambda,lambda_d,maxMuts,s)
#   #If object size
#   while(object.size(tumour)<70){
#     tumour = TumourGrowth(maxsize,beta,lambda,lambda_d,maxMuts, s )
#    }
#   RealVAF = RealDistributionOfVAF(tumour)
#   #save somewhere:
#   save(RealVAF, file = paste("/home/shomeb/p/paalvj/Documents/Masteroppgave/resultater/","RealVAFDist","1mill",beta,lambda,lambda_d,s,9+i,".RData",sep = "_"))
# }
# stopCluster(cl)

  

#######################################################################################








