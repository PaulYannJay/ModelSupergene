#Pkg.add("DataFrames")
#Pkg.add("Requests")
using DataFrames
#using Requests

########### PARAMETERS ---------------------------------------------------------------------

  const NB_GENOTYPES = 128 
  const NB_GENOTYPES_PER_DEME =64 
  const NB_GENOTYPES_PER_CHROM = 8
  const NB_DEMES = 2
  
########### FUNCTIONS ---------------------------------------------------------------------

  ### OBJECTS DESCRIBING THE SYSTEM
  
function createGenotypes()
  
  # genotypeList: characterized by the pairs of alleles on the two chromosomes for the 5^2*2=50 (if 2 demes without dissortative mating) possible genotypes
  # pairs of alleles for each chromosome (indexes 1 and 2): 
  # 1: B-A (inversion)  
  # 2: b-A (inversion)  
  # 3: B-a (inversion)  
  # 4: b-a (inversion)  
  # 5: A-B  
  # 6: A-b
  # 7: a-B
  # 8: a-b
  # (e.g., B-A a-b ---> [1, 8])
  # index 3 : 2 demes: 1 and 2

  genotypeList = ones(Int,NB_GENOTYPES,3) 
  genotypeNb = 1
  for deme in 1:NB_DEMES
    for chrom1 in 1:NB_GENOTYPES_PER_CHROM
        for chrom2 in 1:NB_GENOTYPES_PER_CHROM
       	 genotypeList[genotypeNb,1] = chrom1 # chromosome 1 
       	 genotypeList[genotypeNb,2] = chrom2 # chromosome 2
       	 genotypeList[genotypeNb,3] = deme   # deme
       	 genotypeNb = genotypeNb + 1
        end
    end
  end
  return genotypeList
end

function createPhenotypesEco()
  
  # phenotypeList: phenotypes of the genotypes (cf genotypeList)
  # phenotype for each genotype: 
  # 1: A-B
  # 2: A-b  
  # 3: a-B
  # 4: a-b
  # (e.g., B-A a-b ---> A-B ---> phenotype 1)
  
  phenotypeList = zeros(Int,NB_GENOTYPES) 
  genotypeNb = 1
  for deme in 1:NB_DEMES
    for chrom1 in 1:NB_GENOTYPES_PER_CHROM
      for chrom2 in 1:NB_GENOTYPES_PER_CHROM
        if chrom1==1 || chrom1== 5 || chrom2==1 || chrom2==5 # One of the two allele is BA or AB 
	      phenotypeList[genotypeNb] = 1
        elseif (chrom1==4 || chrom1==8) && (chrom2==4  || chrom2==8) # Both allele are ab or ba
	      phenotypeList[genotypeNb] = 4
        elseif (chrom1==2 || chrom1==6) && chrom2 != 3 && chrom2 != 7 #chrom1 Ab or bA and chrom2 is not Ba or aB --> Ab 
	      phenotypeList[genotypeNb] = 2 # x-x
        elseif (chrom1==1 || chrom1==6) && (chrom2 == 3 || chrom2 == 7) #chrom1 Ab or bA and chrom2 is Ba or aB --> AB 
	      phenotypeList[genotypeNb] = 1 # x-x
        elseif (chrom2==2 || chrom2==6) && chrom1 != 3 && chrom1 != 7 #chrom2 Ab or bA and chrom1 is not Ba or aB --> Ab 
	      phenotypeList[genotypeNb] = 2 # x-x
        elseif (chrom2==2 || chrom2==6) && (chrom1 == 3 || chrom1 == 7) #chrom2 Ab or bA and chrom1 is Ba or aB --> AB 
	      phenotypeList[genotypeNb] = 1 # x-x
		else #Chrom 1 or Chrom 2 is Ba or aB --> Phenotype BA
	      phenotypeList[genotypeNb] = 3 # x-x
        end
        genotypeNb = genotypeNb + 1
      end
    end
  end
  return phenotypeList
end

function createPhenotypesFardeau()
  # phenotypeFardList: phenotypes for Inversion of the genotypes (cf genotypeList)
  # phenotype for each genotype: 
  # 1: Homozygote for the INVERSION 
  # 2: Hetero or without inversion  
  # (e.g., B-A a-b ---> A-B ---> phenotype 1)
  
  phenotypeFardList = zeros(Int,NB_GENOTYPES) 
  genotypeNb = 1
  for deme in 1:NB_DEMES
    for chrom1 in 1:NB_GENOTYPES_PER_CHROM
      for chrom2 in 1:NB_GENOTYPES_PER_CHROM
        if chrom1 < 5 && chrom2 <5
			phenotypeFardList[genotypeNb]=1
		else
			phenotypeFardList[genotypeNb]=2
		end
        genotypeNb = genotypeNb + 1
      end
    end
  end
  return phenotypeFardList
end

function getGenotypeIndexFromGametes()
  
  # genotypeFromGamete: gamete1 + gamete2 -> genotype (index in genotypeList)
  # (e.g., B-A + a-b ---> index of B-A/a-b in genotypeList)
  
  genotypeFromGamete = ones(Int,NB_DEMES,NB_GENOTYPES_PER_CHROM,NB_GENOTYPES_PER_CHROM,1) # 
  genotypeNb = 1
   for deme in 1:NB_DEMES
    for gamete1 in 1:NB_GENOTYPES_PER_CHROM
      for gamete2 in 1:NB_GENOTYPES_PER_CHROM
        genotypeFromGamete[deme, gamete1, gamete2] = genotypeNb # chromosome 1 
        genotypeNb = genotypeNb + 1
      end
    end
  end
  return genotypeFromGamete
end

function createRecombinationMatrix(genotypeList, recombinationRate, recombinationRateRD)

  # formation gametes: recombination rate recombinationRate (except with inversion)
  # chromosome 1 + chromosome 2 --> proportion gamete ( B-A, A-B ... )
  
  proportionGametes = zeros(Float64,NB_GENOTYPES_PER_CHROM,NB_GENOTYPES_PER_CHROM,NB_GENOTYPES_PER_CHROM)  
  for chrom1 in 1:NB_GENOTYPES_PER_CHROM
    for chrom2 in 1:NB_GENOTYPES_PER_CHROM
      if chrom1==chrom2 		# the two chromosome are equals 
        proportionGametes[chrom1,chrom2,chrom1] = 1.0
      elseif (chrom1<5 && chrom2>4) || (chrom2<5 && chrom1>4)	# Individual heterozygote ofr inversion; no recombination
        proportionGametes[chrom1,chrom2,chrom1] = 0.5
        proportionGametes[chrom1,chrom2,chrom2] = 0.5
      elseif (chrom1==5 && chrom2==8) || (chrom1==8 && chrom2==5) # chrom1/chrom2 are A-B/a-b
        proportionGametes[chrom1,chrom2,chrom1] = (1.0-recombinationRate)*0.5
        proportionGametes[chrom1,chrom2,chrom2] = (1.0-recombinationRate)*0.5
        proportionGametes[chrom1,chrom2,6] = recombinationRate*0.5 # A-b recombination
        proportionGametes[chrom1,chrom2,7] = recombinationRate*0.5 # a-B recombination	
      elseif chrom1==5 || chrom1==8 || chrom2==5 || chrom2==8 # launched only if FALSE above; Hom + Het -> recombination is neutral (eg A-B + A-b)
        proportionGametes[chrom1,chrom2,chrom1] = 0.5
        proportionGametes[chrom1,chrom2,chrom2] = 0.5
      elseif (chrom1==6 && chrom2==7) || (chrom1==7 && chrom2==6) # chrom1/chrom2 are A-b/a-B
        proportionGametes[chrom1,chrom2,chrom1] = (1.0-recombinationRate)*0.5
        proportionGametes[chrom1,chrom2,chrom2] = (1.0-recombinationRate)*0.5
        proportionGametes[chrom1,chrom2,5] = recombinationRate*0.5 # A-B recombination
        proportionGametes[chrom1,chrom2,8] = recombinationRate*0.5 # a-b recombination	
      elseif (chrom1==1 && chrom2==4) || (chrom1==4 && chrom2==1) # chrom1/chrom2 are A-B/a-b
        proportionGametes[chrom1,chrom2,chrom1] = (1.0-recombinationRate)*0.5
        proportionGametes[chrom1,chrom2,chrom2] = (1.0-recombinationRate)*0.5
        proportionGametes[chrom1,chrom2,2] = recombinationRate*0.5 # A-b recombination
        proportionGametes[chrom1,chrom2,3] = recombinationRate*0.5 # a-B recombination	
      elseif chrom1==1 || chrom1==4 || chrom2==1 || chrom2==4 # launched only if FALSE above; Hom + Het -> recombination is neutral (eg A-B + A-b)
        proportionGametes[chrom1,chrom2,chrom1] = 0.5
        proportionGametes[chrom1,chrom2,chrom2] = 0.5
      elseif (chrom1==2 && chrom2==3) || (chrom1==3 && chrom2==2) # chrom1/chrom2 are A-b/a-B
        proportionGametes[chrom1,chrom2,chrom1] = (1.0-recombinationRate)*0.5
        proportionGametes[chrom1,chrom2,chrom2] = (1.0-recombinationRate)*0.5
        proportionGametes[chrom1,chrom2,1] = recombinationRate*0.5 # A-B recombination
        proportionGametes[chrom1,chrom2,4] = recombinationRate*0.5 # a-b recombinationend
	  end
    end
  end
  
  # Proportion of each offspring genotype per mating
  
  recombinationMatrix = zeros(Float64,NB_GENOTYPES,NB_GENOTYPES,NB_GENOTYPES) # male + female -> offspring
  getGenotypeFromGametes = getGenotypeIndexFromGametes()
  chrom1=1; chrom2=2
  for mother in 1:NB_GENOTYPES
    propGametesMother = zeros(Float64,NB_GENOTYPES_PER_CHROM) # proportion of gametes produced by the mother
    for gameteMother in 1:NB_GENOTYPES_PER_CHROM
      propGametesMother[gameteMother] = proportionGametes[genotypeList[mother,chrom1],genotypeList[mother,chrom2],gameteMother]
    end
    for father in 1:NB_GENOTYPES 
      if genotypeList[mother,3] == genotypeList[father,3]   #same deme
        deme = genotypeList[mother,3]
        
        propGametesFather = zeros(Float64,NB_GENOTYPES_PER_CHROM) # proportion of gametes produced by the father
        for gameteFather in 1:NB_GENOTYPES_PER_CHROM
          propGametesFather[gameteFather] = proportionGametes[genotypeList[father,chrom1],genotypeList[father,chrom2],gameteFather]
        end
        for gameteMother in 1:NB_GENOTYPES_PER_CHROM
          for gameteFather in 1:NB_GENOTYPES_PER_CHROM
              # resulting offspring, weighted by the probability that the two gametes are produced
              recombinationMatrix[mother,father,getGenotypeFromGametes[deme, gameteMother,gameteFather]] = 
                        recombinationMatrix[mother,father,getGenotypeFromGametes[deme,gameteMother,gameteFather]] + 
                        propGametesFather[gameteFather]*propGametesMother[gameteMother]
          end 
        end        
      end
    end  
  end
  return recombinationMatrix
end

function createSelectionGenotypeVector(selectGoodHomozygotes, selectBadHomozygotes, selectHeterozygotes, additionalSelectInvHomozygotes)
  # selection incurred by each genotype

  phenotypeList=createPhenotypesEco()
  phenotypeFardList=createPhenotypesFardeau()
  selectionGenotypeVector = zeros(Float64,NB_GENOTYPES)
  
  for genotype in 1:NB_GENOTYPES
    deme = genotypeList[genotype,3]
    if phenotypeList[genotype] == 1 	# AB 
      if deme==1
    	if phenotypeFardList[genotype] == 1 	# Inv Hom
        	selectionGenotypeVector[genotype] = selectGoodHomozygotes*additionalSelectInvHomozygotes
		else
        	selectionGenotypeVector[genotype] = selectGoodHomozygotes
		end
      else
    	if phenotypeFardList[genotype] == 1 	# Inv Hom
        	selectionGenotypeVector[genotype] = selectBadHomozygotes*additionalSelectInvHomozygotes
		else
        	selectionGenotypeVector[genotype] = selectBadHomozygotes
		end
      end
    elseif phenotypeList[genotype] == 4 	# ab 
      if deme==1
    	if phenotypeFardList[genotype] == 1 	# Inv Hom
        	selectionGenotypeVector[genotype] = selectBadHomozygotes*additionalSelectInvHomozygotes
		else
        	selectionGenotypeVector[genotype] = selectBadHomozygotes
		end
      else
    	if phenotypeFardList[genotype] == 1 	# Inv Hom
        	selectionGenotypeVector[genotype] = selectGoodHomozygotes*additionalSelectInvHomozygotes
		else
        	selectionGenotypeVector[genotype] = selectGoodHomozygotes
		end
      end 
    elseif phenotypeList[genotype] == 2 || phenotypeList[genotype] == 3 # Ab, aB
    	if phenotypeFardList[genotype] == 1 	# Inv Hom
        	selectionGenotypeVector[genotype] = selectHeterozygotes*additionalSelectInvHomozygotes
		else
        	selectionGenotypeVector[genotype] = selectHeterozygotes
		end
    end
  end
  return selectionGenotypeVector
end

function getChromosomeFrequency(genotypeFrequencies, genotypeList)
  # 1: B-A (inversion)  
  # 2: b-A (inversion)  
  # 3: B-a (inversion)  
  # 4: b-a (inversion)  
  # 5: A-B  
  # 6: A-b
  # 7: a-B
  # 8: a-b
  chromosomeFrequencies = zeros(Float64,NB_GENOTYPES_PER_CHROM)
  for genotype in 1:NB_GENOTYPES
    for chrom in 1:2
      chromosomeFrequencies[genotypeList[genotype,chrom]] = chromosomeFrequencies[genotypeList[genotype,chrom]] +
							    genotypeFrequencies[genotype]
    end
  end
  chromosomeFrequencies = normalizeChromosomeFrequencies(chromosomeFrequencies)
  return chromosomeFrequencies
end

function getChromosomeFrequencyEachDeme(genotypeFrequencies, genotypeList)
  # 1: B-A (inversion)  
  # 2: b-A (inversion)  
  # 3: B-a (inversion)  
  # 4: b-a (inversion)  
  # 5: A-B  
  # 6: A-b
  # 7: a-B
  # 8: a-b
  demeIndex = 3
  chromosomeFrequencies = zeros(Float64,NB_DEMES,NB_GENOTYPES_PER_CHROM)
  for genotype in 1:NB_GENOTYPES
    for chrom in 1:2
      chromosomeFrequencies[genotypeList[genotype,demeIndex], genotypeList[genotype,chrom]] = chromosomeFrequencies[genotypeList[genotype,demeIndex], genotypeList[genotype,chrom]] +
							    genotypeFrequencies[genotype]
    end
  end
  for deme in 1:NB_DEMES
    chromosomeFrequencies[deme,:] = normalizeChromosomeFrequencies(chromosomeFrequencies[deme,:])
  end
  return chromosomeFrequencies
end

function createVectorGenotypeOtherDeme(genotypeList)
  # genotype index corresponding to same genotype in the other deme
  genotypeOtherDemeVector = zeros(Int,NB_GENOTYPES)
  for genotypeSource in 1:NB_GENOTYPES #where individuals are leaving
    found = 0
    genotypeSink=1
    while found==0
        if genotypeList[genotypeSource,1]==genotypeList[genotypeSink,1] && genotypeList[genotypeSource,2]==genotypeList[genotypeSink,2] && genotypeList[genotypeSource,3] != genotypeList[genotypeSink,3]
            #matching genotypes between different demes
            genotypeOtherDemeVector[genotypeSource] = genotypeSink
            found=1
        end    
        genotypeSink=genotypeSink+1
    end
  end
  return genotypeOtherDemeVector
end

function createMutationEffectsMatrix(genotypeList, mutationRate)     # /!\ cumulative probabilities
  mutationEffectsMatrix = zeros(Float64,NB_GENOTYPES,NB_GENOTYPES)
  
  # number of mutations between genotypes from the same deme
  for genotype in 1:NB_GENOTYPES
    for genotypeMutant in 1:NB_GENOTYPES
      if genotypeList[genotype,3]==genotypeList[genotypeMutant,3]  # same deme
	if ((genotypeList[genotype,1]<= 4 && genotypeList[genotypeMutant,1]<= 4)||((genotypeList[genotype,1]>4 && genotypeList[genotypeMutant,1]>4)))&&
	((genotypeList[genotype,2]<= 4 && genotypeList[genotypeMutant,2]<= 4)||((genotypeList[genotype,2]>4 && genotypeList[genotypeMutant,2]>4)))
	  for chrom in 1:2
	    for inv in [0 4] # deal with inverted/non inverted separetely
	    if (genotypeList[genotype,chrom]==1+inv && genotypeList[genotypeMutant,chrom]==2+inv) || 
	       (genotypeList[genotype,chrom]==1+inv && genotypeList[genotypeMutant,chrom]==3+inv) ||
	       (genotypeList[genotype,chrom]==2+inv && genotypeList[genotypeMutant,chrom]==1+inv) ||
	       (genotypeList[genotype,chrom]==2+inv && genotypeList[genotypeMutant,chrom]==4+inv) ||
	       (genotypeList[genotype,chrom]==3+inv && genotypeList[genotypeMutant,chrom]==1+inv) ||
	       (genotypeList[genotype,chrom]==3+inv && genotypeList[genotypeMutant,chrom]==4+inv) ||
	       (genotypeList[genotype,chrom]==4+inv && genotypeList[genotypeMutant,chrom]==2+inv) ||
	       (genotypeList[genotype,chrom]==4+inv && genotypeList[genotypeMutant,chrom]==3+inv)
	      mutationEffectsMatrix[genotype,genotypeMutant]=mutationEffectsMatrix[genotype,genotypeMutant]+1 #one mutation
	    elseif (genotypeList[genotype,chrom]==1+inv && genotypeList[genotypeMutant,chrom]==4+inv) || 
	       (genotypeList[genotype,chrom]==2+inv && genotypeList[genotypeMutant,chrom]==3+inv) ||
	       (genotypeList[genotype,chrom]==3+inv && genotypeList[genotypeMutant,chrom]==2+inv) ||
	       (genotypeList[genotype,chrom]==4+inv && genotypeList[genotypeMutant,chrom]==1+inv)
	      mutationEffectsMatrix[genotype,genotypeMutant]=mutationEffectsMatrix[genotype,genotypeMutant]+2 #two mutations
	    end
	    end
	  end
	end
      if mutationEffectsMatrix[genotype,genotypeMutant]>0.0 && mutationRate>0.0
	mutationEffectsMatrix[genotype,genotypeMutant] = mutationRate^mutationEffectsMatrix[genotype,genotypeMutant]
      else
        mutationEffectsMatrix[genotype,genotypeMutant] = 0.0
      end
      end
    end
#        println(sum(mutationEffectsMatrix[genotype,1:end]))
  end
#   println(mutationEffectsMatrix[37,1:end])
  # proba to mutate to other genotypes (1/(number of mutations)); mutation rate per individual cf. other function
  for genotype in 1:NB_GENOTYPES
    sumOneGenotype = sum(mutationEffectsMatrix[genotype,1:end])
    mutationEffectsMatrix[genotype,1:end] = mutationEffectsMatrix[genotype,1:end]/sumOneGenotype
    # cumulative probabilities
    for newMutant in 2:NB_GENOTYPES
      mutationEffectsMatrix[genotype,newMutant] = mutationEffectsMatrix[genotype,newMutant-1] + mutationEffectsMatrix[genotype,newMutant]
    end    
  end
          
#   println(mutationEffectsMatrix[2,1:NB_GENOTYPES])
  return mutationEffectsMatrix
end


function calculateMutationRatePerIndividual(genotypeList, mutationRate)     # ~same function than above
  mutationEffectsMatrix = zeros(Float64,NB_GENOTYPES,NB_GENOTYPES)
  
  # number of mutations between genotypes from the same deme
  for genotype in 1:NB_GENOTYPES
    for genotypeMutant in 1:NB_GENOTYPES
      if genotypeList[genotype,3]==genotypeList[genotypeMutant,3]  # same deme
	if ((genotypeList[genotype,1]<= 4 && genotypeList[genotypeMutant,1]<= 4)||((genotypeList[genotype,1]>4 && genotypeList[genotypeMutant,1]>4)))&&
	((genotypeList[genotype,2]<= 4 && genotypeList[genotypeMutant,2]<= 4)||((genotypeList[genotype,2]>4 && genotypeList[genotypeMutant,2]>4)))
	  for chrom in 1:2
	    for inv in [0 4] # deal with inverted/non inverted separetely
	    if (genotypeList[genotype,chrom]==1+inv && genotypeList[genotypeMutant,chrom]==2+inv) || 
	       (genotypeList[genotype,chrom]==1+inv && genotypeList[genotypeMutant,chrom]==3+inv) ||
	       (genotypeList[genotype,chrom]==2+inv && genotypeList[genotypeMutant,chrom]==1+inv) ||
	       (genotypeList[genotype,chrom]==2+inv && genotypeList[genotypeMutant,chrom]==4+inv) ||
	       (genotypeList[genotype,chrom]==3+inv && genotypeList[genotypeMutant,chrom]==1+inv) ||
	       (genotypeList[genotype,chrom]==3+inv && genotypeList[genotypeMutant,chrom]==4+inv) ||
	       (genotypeList[genotype,chrom]==4+inv && genotypeList[genotypeMutant,chrom]==2+inv) ||
	       (genotypeList[genotype,chrom]==4+inv && genotypeList[genotypeMutant,chrom]==3+inv)
	      mutationEffectsMatrix[genotype,genotypeMutant]=mutationEffectsMatrix[genotype,genotypeMutant]+1 #one mutation
	    elseif (genotypeList[genotype,chrom]==1+inv && genotypeList[genotypeMutant,chrom]==4+inv) || 
	       (genotypeList[genotype,chrom]==2+inv && genotypeList[genotypeMutant,chrom]==3+inv) ||
	       (genotypeList[genotype,chrom]==3+inv && genotypeList[genotypeMutant,chrom]==2+inv) ||
	       (genotypeList[genotype,chrom]==4+inv && genotypeList[genotypeMutant,chrom]==1+inv)
	      mutationEffectsMatrix[genotype,genotypeMutant]=mutationEffectsMatrix[genotype,genotypeMutant]+2 #two mutations
	    end
	    end
	  end
	end
      if mutationEffectsMatrix[genotype,genotypeMutant]>0.0 && mutationRate>0.0
	mutationEffectsMatrix[genotype,genotypeMutant] = mutationRate^mutationEffectsMatrix[genotype,genotypeMutant]
      else
        mutationEffectsMatrix[genotype,genotypeMutant] = 0.0
      end
      end
    end
  end
  mutationRatePerIndividual = sum(mutationEffectsMatrix[1,1:end])
  return mutationRatePerIndividual
end


function createInversionEffectsMatrixTwoWays(genotypeList, invasionRate)     # /!\ cumulative probabilities
  inversionEffectsMatrix = zeros(Float64,NB_GENOTYPES,NB_GENOTYPES)
  
  # number of inversions between genotypes from the same deme
  for genotype in 1:NB_GENOTYPES
    for genotypeMutant in 1:NB_GENOTYPES
      if genotypeList[genotype,3]==genotypeList[genotypeMutant,3]  # same deme
            if max(genotypeList[genotype,1],genotypeList[genotypeMutant,1])-min(genotypeList[genotype,1],genotypeList[genotypeMutant,1])==4
                if max(genotypeList[genotype,2],genotypeList[genotypeMutant,2])-min(genotypeList[genotype,2],genotypeList[genotypeMutant,2])==4
                    inversionEffectsMatrix[genotype,genotypeMutant]=inversionEffectsMatrix[genotype,genotypeMutant]+2 #two inversion
                elseif genotypeList[genotype,2]==genotypeList[genotypeMutant,2]
                    inversionEffectsMatrix[genotype,genotypeMutant]=inversionEffectsMatrix[genotype,genotypeMutant]+1 #one inversion
                end
            elseif genotypeList[genotype,1]==genotypeList[genotypeMutant,1]
                if max(genotypeList[genotype,2],genotypeList[genotypeMutant,2])-min(genotypeList[genotype,2],genotypeList[genotypeMutant,2])==4
                    inversionEffectsMatrix[genotype,genotypeMutant]=inversionEffectsMatrix[genotype,genotypeMutant]+1 #one inversion
                end
            end
      if inversionEffectsMatrix[genotype,genotypeMutant]>0 && invasionRate>0
        inversionEffectsMatrix[genotype,genotypeMutant] = invasionRate^inversionEffectsMatrix[genotype,genotypeMutant]
      else
        inversionEffectsMatrix[genotype,genotypeMutant] = 0.0
      end              
      end
    end
#        println(sum(inversionEffectsMatrix[genotype,1:end]))
  end
#   println(inversionEffectsMatrix[1,1:end])
  # proba to mutate to other genotypes (1/(number of mutations)); mutation rate per individual cf. other function
  for genotype in 1:NB_GENOTYPES
    sumOneGenotype = sum(inversionEffectsMatrix[genotype,1:end])
    inversionEffectsMatrix[genotype,1:end] = inversionEffectsMatrix[genotype,1:end]/sumOneGenotype
    # cumulative probabilities
    for newMutant in 2:NB_GENOTYPES
      inversionEffectsMatrix[genotype,newMutant] = inversionEffectsMatrix[genotype,newMutant-1] + inversionEffectsMatrix[genotype,newMutant]
    end    
  end
# #   println(inversionEffectsMatrix[2,1:NB_GENOTYPES])
  return inversionEffectsMatrix
end

function calculateInversionRatePerIndividualTwoWays(genotypeList, invasionRate)     # ~same function than above
  inversionEffectsMatrix = zeros(Float64,NB_GENOTYPES,NB_GENOTYPES)
  
  # number of inversions between genotypes from the same deme
  for genotype in 1:NB_GENOTYPES
    for genotypeMutant in 1:NB_GENOTYPES
      if genotypeList[genotype,3]==genotypeList[genotypeMutant,3]  # same deme
            if max(genotypeList[genotype,1],genotypeList[genotypeMutant,1])-min(genotypeList[genotype,1],genotypeList[genotypeMutant,1])==4
                if max(genotypeList[genotype,2],genotypeList[genotypeMutant,2])-min(genotypeList[genotype,2],genotypeList[genotypeMutant,2])==4
                    inversionEffectsMatrix[genotype,genotypeMutant]=inversionEffectsMatrix[genotype,genotypeMutant]+2.0 #two inversion
                elseif genotypeList[genotype,2]==genotypeList[genotypeMutant,2]
                    inversionEffectsMatrix[genotype,genotypeMutant]=inversionEffectsMatrix[genotype,genotypeMutant]+1.0 #one inversion
                end
            elseif genotypeList[genotype,1]==genotypeList[genotypeMutant,1]
                if max(genotypeList[genotype,2],genotypeList[genotypeMutant,2])-min(genotypeList[genotype,2],genotypeList[genotypeMutant,2])==4
                    inversionEffectsMatrix[genotype,genotypeMutant]=inversionEffectsMatrix[genotype,genotypeMutant]+1.0 #one inversion
                end
            end
      end
      if inversionEffectsMatrix[genotype,genotypeMutant]>0 && invasionRate>0
        inversionEffectsMatrix[genotype,genotypeMutant] = invasionRate^inversionEffectsMatrix[genotype,genotypeMutant]
      else
        inversionEffectsMatrix[genotype,genotypeMutant] = 0.0
      end   
    end 
#        println(sum(inversionEffectsMatrix[genotype,1:end]))
  end
  inversionRatePerIndividual = zeros(Float64,NB_GENOTYPES)
  for genotype in 1:NB_GENOTYPES
    inversionRatePerIndividual[genotype] = sum(inversionEffectsMatrix[genotype,1:end])
  end
  return inversionRatePerIndividual
end


function createInversionEffectsMatrixOneWay(genotypeList, invasionRate)     # /!\ cumulative probabilities  INVERSION IN ONE WAY ONLY
  inversionEffectsMatrix = zeros(Float64,NB_GENOTYPES,NB_GENOTYPES)
  
  # number of inversions between genotypes from the same deme
  for genotype in 1:NB_GENOTYPES
    for genotypeMutant in 1:NB_GENOTYPES
      if genotypeList[genotype,3]==genotypeList[genotypeMutant,3]  # same deme
            if genotypeList[genotype,1]-genotypeList[genotypeMutant,1]==4
                if genotypeList[genotype,2]-genotypeList[genotypeMutant,2]==4
                    inversionEffectsMatrix[genotype,genotypeMutant]=inversionEffectsMatrix[genotype,genotypeMutant]+2.0 #two inversion
                elseif genotypeList[genotype,2]==genotypeList[genotypeMutant,2]
                    inversionEffectsMatrix[genotype,genotypeMutant]=inversionEffectsMatrix[genotype,genotypeMutant]+1.0 #one inversion
                end
            elseif genotypeList[genotype,1]==genotypeList[genotypeMutant,1]
                if genotypeList[genotype,2]-genotypeList[genotypeMutant,2]==4
                    inversionEffectsMatrix[genotype,genotypeMutant]=inversionEffectsMatrix[genotype,genotypeMutant]+1.0 #one inversion
                end
            end
      if inversionEffectsMatrix[genotype,genotypeMutant]>0 && invasionRate>0
        inversionEffectsMatrix[genotype,genotypeMutant] = invasionRate^inversionEffectsMatrix[genotype,genotypeMutant]
      else
        inversionEffectsMatrix[genotype,genotypeMutant] = 0.0
      end                      
      end
    end
#        println(sum(inversionEffectsMatrix[genotype,1:end]))
  end
#   println(inversionEffectsMatrix[1,1:end])
  # proba to mutate to other genotypes (1/(number of mutations)); mutation rate per individual cf. other function
  for genotype in 1:NB_GENOTYPES
    sumOneGenotype = sum(inversionEffectsMatrix[genotype,1:end])
    inversionEffectsMatrix[genotype,1:end] = inversionEffectsMatrix[genotype,1:end]/sumOneGenotype
    # cumulative probabilities
    for newMutant in 2:NB_GENOTYPES
      inversionEffectsMatrix[genotype,newMutant] = inversionEffectsMatrix[genotype,newMutant-1] + inversionEffectsMatrix[genotype,newMutant]
    end    
  end
# #   println(inversionEffectsMatrix[2,1:NB_GENOTYPES])
  return inversionEffectsMatrix
end

function calculateInversionRatePerIndividualOneWay(genotypeList, invasionRate)     # ~same function than above INVERSION IN ONE WAY ONLY
  inversionEffectsMatrix = zeros(Float64,NB_GENOTYPES,NB_GENOTYPES)
  
  # number of inversions between genotypes from the same deme
  for genotype in 1:NB_GENOTYPES
    for genotypeMutant in 1:NB_GENOTYPES
      if genotypeList[genotype,3]==genotypeList[genotypeMutant,3]  # same deme
            if genotypeList[genotype,1]-genotypeList[genotypeMutant,1]==4
                if genotypeList[genotype,2]-genotypeList[genotypeMutant,2]==4
                    inversionEffectsMatrix[genotype,genotypeMutant]=inversionEffectsMatrix[genotype,genotypeMutant]+2.0 #two inversion
                elseif genotypeList[genotype,2]==genotypeList[genotypeMutant,2]
                    inversionEffectsMatrix[genotype,genotypeMutant]=inversionEffectsMatrix[genotype,genotypeMutant]+1.0 #one inversion
                end
            elseif genotypeList[genotype,1]==genotypeList[genotypeMutant,1]
                if genotypeList[genotype,2]-genotypeList[genotypeMutant,2]==4
                    inversionEffectsMatrix[genotype,genotypeMutant]=inversionEffectsMatrix[genotype,genotypeMutant]+1.0 #one inversion
                end
            end
      end
      if inversionEffectsMatrix[genotype,genotypeMutant]>0 && invasionRate>0
        inversionEffectsMatrix[genotype,genotypeMutant] = invasionRate^inversionEffectsMatrix[genotype,genotypeMutant]
      else
        inversionEffectsMatrix[genotype,genotypeMutant] = 0.0
      end   
    end 
#        println(sum(inversionEffectsMatrix[genotype,1:end]))
  end

  inversionRatePerIndividual = zeros(Float64,NB_GENOTYPES)
  for genotype in 1:NB_GENOTYPES
    inversionRatePerIndividual[genotype] = sum(inversionEffectsMatrix[genotype,1:end])
  end

  return inversionRatePerIndividual
end

  ### INITIALIZATION ## Different initial state possible

function initializeGenotypeFrequencies(genotypeList, typeOfInitialization, percentBAIni, recombinationMatrix, selectionGenotypeVector, migrationRate, migrationUnidirectionnal)
  if typeOfInitialization=="Uniform"
  
    genotypeFrequencies = 1/NB_GENOTYPES_PER_DEME * ones(Float64,NB_GENOTYPES)
elseif typeOfInitialization=="OnlyBA"
  
    #OnlyAB initially
    genotypeFrequencies = zeros(Float64,NB_GENOTYPES)
    sumGenotypesOnlyAB = 0.0
    for i in 1:NB_GENOTYPES
      if genotypeList[i,1] == 1 && genotypeList[i,2] == 1  # Hom B-A
        sumGenotypesOnlyAB=sumGenotypesOnlyAB+1.0
      end
    end    
    for i in 1:NB_GENOTYPES
      if genotypeList[i,1] == 1 && genotypeList[i,2] == 1  # Hom B-A
        genotypeFrequencies[i] = 1/sumGenotypesOnlyAB
      end
    end
    genotypeFrequencies = normalizeGenotypeFrequenciesTwoDemes(genotypeFrequencies)
  elseif typeOfInitialization=="Onlyab"
  
    #OnlyAB initially
    genotypeFrequencies = zeros(Float64,NB_GENOTYPES)
    sumGenotypesOnlyAB = 0.0
    for i in 1:NB_GENOTYPES
      if genotypeList[i,1] == 8 && genotypeList[i,2] == 8  # Hom A-B
        sumGenotypesOnlyAB=sumGenotypesOnlyAB+1.0
      end
    end    
    for i in 1:NB_GENOTYPES
      if genotypeList[i,1] == 8 && genotypeList[i,2] == 8  # Hom A-B
        genotypeFrequencies[i] = 1/sumGenotypesOnlyAB
      end
    end
    genotypeFrequencies = normalizeGenotypeFrequenciesTwoDemes(genotypeFrequencies)
 
  elseif typeOfInitialization=="OnlyAB"
  
    #OnlyAB initially
    genotypeFrequencies = zeros(Float64,NB_GENOTYPES)
    sumGenotypesOnlyAB = 0.0
    for i in 1:NB_GENOTYPES
      if genotypeList[i,1] == 5 && genotypeList[i,2] == 5  # Hom A-B
        sumGenotypesOnlyAB=sumGenotypesOnlyAB+1.0
      end
    end    
    for i in 1:NB_GENOTYPES
      if genotypeList[i,1] == 5 && genotypeList[i,2] == 5  # Hom A-B
        genotypeFrequencies[i] = 1/sumGenotypesOnlyAB
      end
    end
    genotypeFrequencies = normalizeGenotypeFrequenciesTwoDemes(genotypeFrequencies)
  elseif typeOfInitialization=="OnlyABab"
  
    #Only AB and ab initially
    genotypeFrequencies = zeros(Float64,NB_GENOTYPES)
    sumGenotypesOnlyABab = 0.0
    for i in 1:NB_GENOTYPES
      if genotypeList[i,1] == 5 && genotypeList[i,2] == 5  # Hom A-B
        sumGenotypesOnlyABab=sumGenotypesOnlyABab+1.0
      elseif genotypeList[i,1] == 8 && genotypeList[i,2] == 8  # Hom a-b
        sumGenotypesOnlyABab=sumGenotypesOnlyABab+1.0
      end
    end    
    for i in 1:NB_GENOTYPES
      if genotypeList[i,1] == 5 && genotypeList[i,2] == 5  # Hom A-B
        genotypeFrequencies[i] = 1/sumGenotypesOnlyABab
      elseif genotypeList[i,1] == 8 && genotypeList[i,2] == 8  # Hom a-b
        genotypeFrequencies[i] = 1/sumGenotypesOnlyABab
      end
    end
    genotypeFrequencies = normalizeGenotypeFrequenciesTwoDemes(genotypeFrequencies)
  elseif typeOfInitialization=="UniformNoBA"
  
    #no BA initially
    genotypeFrequencies = zeros(Float64,NB_GENOTYPES)
    sumGenotypesNoBA = 0.0
    for i in 1:NB_GENOTYPES
      if genotypeList[i,1] >4 && genotypeList[i,2] > 4  # inversion B-A absent
        sumGenotypesNoBA = sumGenotypesNoBA + 1.0
      end
    end
    for i in 1:NB_GENOTYPES
      if genotypeList[i,1] >4 && genotypeList[i,2] >4  # inversion B-A absent
        genotypeFrequencies[i] = 1/sumGenotypesNoBA
      end
    end
    genotypeFrequencies = normalizeGenotypeFrequenciesTwoDemes(genotypeFrequencies)
    
  elseif typeOfInitialization=="InvasionBA"
  
    #OnlyAB initially
    genotypeFrequencies = zeros(Float64,NB_GENOTYPES)
    sumGenotypesOnlyAB = 0.0
    for i in 1:NB_GENOTYPES
      if genotypeList[i,1] == 5 && genotypeList[i,2] == 5  # Hom A-B
        sumGenotypesOnlyAB=sumGenotypesOnlyAB+1.0
      end
    end    
    for i in 1:NB_GENOTYPES
      if genotypeList[i,1] == 5 && genotypeList[i,2] == 5  # Hom A-B
        genotypeFrequencies[i] = 1/sumGenotypesOnlyAB
      end
    end
    genotypeFrequencies = normalizeGenotypeFrequenciesTwoDemes(genotypeFrequencies)

      #equilibrium reached after Tsim generations
    TsimBefore = 200
    genotypeFrequencies = freqAfterSimulation(genotypeList, genotypeFrequencies, recombinationMatrix, selectionGenotypeVector, migrationRate, migrationUnidirectionnal, TsimBefore)

      #implementation of BA from AB
    genotypeFrequencies = implementBAgenotype(genotypeList, genotypeFrequencies, percentBAIni)
    
    genotypeFrequencies = normalizeGenotypeFrequenciesTwoDemes(genotypeFrequencies)

  end
  return genotypeFrequencies
end

function implementBAgenotype(genotypeList, genotypeFrequencies, percentBAIni)
 
    #implementation of BA from AB
  alleleBA = 1
  alleleAB = 5
  chrom1=1
  chrom2=2
  demeIndex = 3
  for i in 1:NB_GENOTYPES
    if genotypeList[i,chrom1] == alleleAB && genotypeList[i,chrom2] == alleleAB  # A-B on the two chromosomes
        for j in 1:NB_GENOTYPES
            if genotypeList[j,demeIndex] == genotypeList[i,demeIndex] && genotypeList[j,chrom1] == alleleBA && genotypeList[j,chrom2] == alleleBA 
                genotypeFrequencies[j] = percentBAIni * genotypeFrequencies[i] #B-A Hom
                genotypeFrequencies[i] = (1.0-percentBAIni) * genotypeFrequencies[i]
      
            end
        end
    elseif genotypeList[i,chrom1] == alleleAB # A-B on the first chromosome
      for j in 1:NB_GENOTYPES
        if genotypeList[j,demeIndex] == genotypeList[i,demeIndex] && genotypeList[j,chrom2] == genotypeList[i,chrom2] && genotypeList[j,chrom1] == alleleBA 
            # same but with B-A on the first chromosome
            genotypeFrequencies[j] = percentBAIni * genotypeFrequencies[i]
            genotypeFrequencies[i] = (1.0-percentBAIni) * genotypeFrequencies[i]
        end
      end
    elseif genotypeList[i,chrom2] == alleleAB # A-B on the second chromosome
      for j in 1:NB_GENOTYPES
        if genotypeList[j,demeIndex] == genotypeList[i,demeIndex] && genotypeList[j,chrom1] == genotypeList[i,chrom1] && genotypeList[j,chrom2] == alleleBA 
            # same but with B-A on the second chromosome
            genotypeFrequencies[j] = percentBAIni * genotypeFrequencies[i]
            genotypeFrequencies[i] = (1.0-percentBAIni) * genotypeFrequencies[i]
        end
      end
    end
  end
  
  genotypeFrequencies = normalizeGenotypeFrequenciesTwoDemes(genotypeFrequencies)
  return genotypeFrequencies
end

  ### MIGRATION

function getFrequenciesAfterMigrationParapatry(genotypeList, genotypeFrequencies, genotypeOtherDemeVector, migrationRate, migrationUnidirectionnal)

  # genotype frequencies after migration has occurred 
  frequenciesAfterMigration = zeros(Float64,NB_GENOTYPES)

  if migrationUnidirectionnal==false
    for genotypeSource in 1:NB_GENOTYPES #where individuals are leaving
        genotypeSink = genotypeOtherDemeVector[genotypeSource]
        frequenciesAfterMigration[genotypeSource] = frequenciesAfterMigration[genotypeSource] + (1.0-migrationRate) * genotypeFrequencies[genotypeSource]
        frequenciesAfterMigration[genotypeSink] =  frequenciesAfterMigration[genotypeSink] + migrationRate * genotypeFrequencies[genotypeSource]
    end
    frequenciesAfterMigration = normalizeGenotypeFrequenciesTwoDemes(frequenciesAfterMigration)
  else
    for genotypeSource in 1:NB_GENOTYPES #where individuals are leaving
        if genotypeList[genotypeSource,3]==2
        genotypeSink = genotypeOtherDemeVector[genotypeSource]
        frequenciesAfterMigration[genotypeSource] = frequenciesAfterMigration[genotypeSource] + (1.0-migrationRate) * genotypeFrequencies[genotypeSource]
        frequenciesAfterMigration[genotypeSink] =  frequenciesAfterMigration[genotypeSink] + migrationRate * genotypeFrequencies[genotypeSource]
        else
        genotypeSink = genotypeOtherDemeVector[genotypeSource]
        frequenciesAfterMigration[genotypeSource] = frequenciesAfterMigration[genotypeSource] + genotypeFrequencies[genotypeSource]
        end
    end
    frequenciesAfterMigration = normalizeGenotypeFrequenciesTwoDemes(frequenciesAfterMigration)
  end

  return frequenciesAfterMigration
end

  ### REPRODUCTION

function getMatingFraction(genotypeFrequencies, disassortativeLevel)

  # fraction of mating between each genotype
  # here random mating, but we can implement dissortative mating later
  
  fractionMating = zeros(Float64,NB_GENOTYPES,NB_GENOTYPES) # row: female; column: male

  # overall fraction of mating between a given female and a given male
  for ind1 in 1:NB_GENOTYPES
    for ind2 in ind1:NB_GENOTYPES
      if genotypeList[ind1,3] == genotypeList[ind2,3] # same deme;mating do not occur if individuals from different deme!
        fractionMating[ind1,ind2] = genotypeFrequencies[ind1] * genotypeFrequencies[ind2] # ind1 = mother
        fractionMating[ind2,ind1] = genotypeFrequencies[ind2] * genotypeFrequencies[ind1] # ind2 = mother
      end
    end
  end

  return fractionMating
end

function getFrequenciesOffspringDeterministic(genotypeList, genotypeFrequencies, fractionMating, recombinationMatrix)
  
  # genotype frequencies in offspring
  genotypeFrequenciesOffspring = zeros(Float64,NB_GENOTYPES)
  
  for ind2 in 1:NB_GENOTYPES
    for ind1 in ind2:NB_GENOTYPES
      if genotypeList[ind1,3] == genotypeList[ind2,3] # same deme
        for offspring in 1:NB_GENOTYPES
            if genotypeList[offspring,3] == genotypeList[ind1,3] # same deme
                if ind1!=ind2
                    genotypeFrequenciesOffspring[offspring] = genotypeFrequenciesOffspring[offspring] + 
                                            2.0*fractionMating[ind1,ind2]*recombinationMatrix[ind1,ind2,offspring]     
                else
                    genotypeFrequenciesOffspring[offspring] = genotypeFrequenciesOffspring[offspring] + 
                                            fractionMating[ind1,ind2]*recombinationMatrix[ind1,ind2,offspring]               
                end
            end
        end
      end
    end
  end
  genotypeFrequenciesOffspring = normalizeGenotypeFrequenciesTwoDemes(genotypeFrequenciesOffspring)

  return genotypeFrequenciesOffspring
end

  ### VIABILITY SELECTION 
  
function getFrequenciesAfterSelection(genotypeList, genotypeFrequencies, selectionGenotypeVector)

  # genotype frequencies after viability selection has occured
  frequenciesAfterSelection = zeros(Float64,NB_GENOTYPES)
  
  for genotype in 1:NB_GENOTYPES
    frequenciesAfterSelection[genotype] = genotypeFrequencies[genotype] * selectionGenotypeVector[genotype]
  end
  frequenciesAfterSelection = normalizeGenotypeFrequenciesTwoDemes(frequenciesAfterSelection)
  
  return frequenciesAfterSelection
end

  ### NORMALIZATION

function normalizeChromosomeFrequencies(chromosomeFrequencies)
  # normalize to ensure sum(chromosomeFrequencies) = 1  (like old version of normalizeGenotypeFrequencies, cf above)
  genotypeChromosomeNormalized = chromosomeFrequencies / sum(chromosomeFrequencies)
  return genotypeChromosomeNormalized
end

function normalizeGenotypeFrequenciesTwoDemes(genotypeFrequencies)
  # normalize to ensure sum(genotypeFrequencies) = 1 in each deme (!)
  sumFreqEachDeme = zeros(Float64,2)
  genotypeFrequenciesNormalized = zeros(Float64,NB_GENOTYPES)
  for genotype in 1:NB_GENOTYPES
    sumFreqEachDeme[genotypeList[genotype,3]] = sumFreqEachDeme[genotypeList[genotype,3]] + genotypeFrequencies[genotype]
  end
  for genotype in 1:NB_GENOTYPES
    genotypeFrequenciesNormalized[genotype] = genotypeFrequencies[genotype]/sumFreqEachDeme[genotypeList[genotype,3]]
  end
  return genotypeFrequenciesNormalized
end

function normalizeGenotypeFrequenciesTotal(genotypeFrequencies)
  # normalize to ensure sum(genotypeFrequencies) = 1 in total (no distinction between demes!)
  sumFreq = 0.0
  genotypeFrequenciesNormalized = zeros(Float64,NB_GENOTYPES)
  for genotype in 1:NB_GENOTYPES
    sumFreq = sumFreq + genotypeFrequencies[genotype]
  end
  for genotype in 1:NB_GENOTYPES
    genotypeFrequenciesNormalized[genotype] = genotypeFrequencies[genotype]/sumFreq
  end
  return genotypeFrequenciesNormalized
end

  ### GENETIC DRIFT WITHIN FINITE POPULATION

function getFrequenciesFinitePopulation(genotypeFrequencies, carryingCapacity, mutationEffectsMatrix, mutationRatePerIndividual, inversionEffectsMatrix, inversionRatePerIndividual)
  # We sample individuals from frequency distribution (no distinction between demes)
  # i.e., total carrying capacity
  genotypeFrequencies = normalizeGenotypeFrequenciesTotal(genotypeFrequencies)
  genotypeFrequenciesCumulative = zeros(Float64,NB_GENOTYPES)
  genotypeFrequenciesCumulative[1] = genotypeFrequencies[1]
  for genotype in 2:NB_GENOTYPES
    genotypeFrequenciesCumulative[genotype] = genotypeFrequenciesCumulative[genotype-1] + genotypeFrequencies[genotype]
    if genotypeFrequenciesCumulative[genotype]>1.0
      genotypeFrequenciesCumulative[genotype] = 1.0
    end
  end
  
  genotypeFrequenciesRandom = zeros(Float64,NB_GENOTYPES)
  for ind in 1:carryingCapacity
    randNb = rand()
    genotype = 1
    Found = false
    while genotype<=NB_GENOTYPES && Found == false
      if randNb<genotypeFrequenciesCumulative[genotype]
        Found = true
        
        genotypeFinal = copy(genotype)

        if rand()<mutationRatePerIndividual  # at least one mutation occurs
            randNb2 = rand() # what mutation?
            genotype2 = 1
            Found2 = false 
            while genotype2<NB_GENOTYPES && Found2 == false
                if randNb2<mutationEffectsMatrix[genotypeFinal,genotype2]
                Found2 = true
                genotypeFinal = copy(genotype2)
                end
                genotype2 = genotype2+1
            end
        end

        if rand()<inversionRatePerIndividual[genotypeFinal]  # at least one inversion occurs
            randNb3 = rand() # what inversion?
            genotype3 = 1
            Found3 = false
            while genotype3<NB_GENOTYPES && Found3 == false
                if randNb3<inversionEffectsMatrix[genotypeFinal,genotype3]
                Found3 = true
                genotypeFinal = copy(genotype3)
                end
                genotype3 = genotype3+1
            end
        end
        genotypeFrequenciesRandom[genotypeFinal] = genotypeFrequenciesRandom[genotypeFinal] + 1.0
        
      end
      genotype = genotype+1 
    end
  end
  genotypeFrequenciesRandom = normalizeGenotypeFrequenciesTwoDemes(genotypeFrequenciesRandom)
  return genotypeFrequenciesRandom
end

    ## Remove gentype At very Low Frequencies

function RemoveLowFrequencies(genotypeFrequencies)
	
  frequenciesAfterPurge = zeros(Float64,NB_GENOTYPES)
  for genotype in 1:NB_GENOTYPES
  	if genotypeFrequencies[genotype] < 0.0001
  		frequenciesAfterPurge[genotype]=0.0
	else
		frequenciesAfterPurge[genotype]=genotypeFrequencies[genotype]
	end
  end
  frequenciesAfterPurge=normalizeGenotypeFrequenciesTwoDemes(frequenciesAfterPurge)
  return frequenciesAfterPurge
end

  ### SIMULATION
  
function freqAfterSimulation(genotypeList, genotypeFrequencies, recombinationMatrix, mutationEffectsMatrix, mutationRatePerIndividual, inversionEffectsMatrix, inversionRatePerIndividual, selectionGenotypeVector, migrationRate, migrationUnidirectionnal, Tsim)

  for i in 1:Tsim
    # migration
    genotypeFrequencies = getFrequenciesAfterMigrationParapatry(genotypeList, genotypeFrequencies, genotypeOtherDemeVector, migrationRate, migrationUnidirectionnal)
    
    # reproduction
    fractionMating = getMatingFraction(genotypeFrequencies, disassortativeLevel)
    genotypeFrequencies = getFrequenciesOffspringDeterministic(genotypeList, genotypeFrequencies, fractionMating, recombinationMatrix)
    
    # selection
    genotypeFrequencies = getFrequenciesAfterSelection(genotypeList, genotypeFrequencies, selectionGenotypeVector)

    # genetic drift
    if boolFinitePopulation==1
      genotypeFrequencies = getFrequenciesFinitePopulation(genotypeFrequencies, carryingCapacity, mutationEffectsMatrix, mutationRatePerIndividual, inversionEffectsMatrix, inversionRatePerIndividual)
    end
    
  end

  return genotypeFrequencies
end

########### ---------------------------------------------------------------------
