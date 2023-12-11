using DataFrames
#using Requests

########### OUTPUT ---------------------------------------------------------------------


  ### PLOT CSV
  
function csvOneSimulation(genotypeList, recombinationMatrix, disassortativeMatingScript, disassortativeLevel, selectionGenotypeVector, genotypeOtherDemeVector, migrationRate, migrationUnidirectionnal, boolFinitePopulation, carryingCapacity, mutationRate, inversionRate, inversionType, TstopMigration,Tend, nameFile)

  if isfile("Data/SingleTimeSeries"*nameFile*".txt")==true
    rm("Data/SingleTimeSeries"*nameFile*".txt")
  end
  file = open("Data/SingleTimeSeries"*nameFile*".txt", "w")
  write(file,"time;genotype;frequency \n")
  if isfile("Data/SingleTimeSeries2"*nameFile*".txt")==true
    rm("Data/SingleTimeSeries2"*nameFile*".txt")
  end
  file2 = open("Data/SingleTimeSeries2"*nameFile*".txt", "w")
  write(file2,"time;deme;chromosome;frequency \n")
  
  println(nameFile)
  
  #mutation effects
  mutationEffectsMatrix = createMutationEffectsMatrix(genotypeList, mutationRate)
  mutationRatePerIndividual = calculateMutationRatePerIndividual(genotypeList, mutationRate) 

  #inversion effects
  if inversionType==0
    inversionEffectsMatrix = createInversionEffectsMatrixTwoWays(genotypeList, inversionRate)
    inversionRatePerIndividual = calculateInversionRatePerIndividualTwoWays(genotypeList, inversionRate) 
  elseif inversionType==1
    inversionEffectsMatrix = createInversionEffectsMatrixOneWay(genotypeList, inversionRate)
    inversionRatePerIndividual = calculateInversionRatePerIndividualOneWay(genotypeList, inversionRate) 
  end
  
  genotypeFrequencies = initializeGenotypeFrequencies(genotypeList, typeOfInitialization, percentBAIni, recombinationMatrix, selectionGenotypeVector, migrationRate, migrationUnidirectionnal)
    
  t = 0
  while t<=Tend
    for genotype in 1:NB_GENOTYPES
      write(file, "$(t);$(genotype);$(genotypeFrequencies[genotype]) \n")
    end
    chromosomeFrequencies = getChromosomeFrequencyEachDeme(genotypeFrequencies, genotypeList)
    for deme in 1:NB_DEMES
        for chromType in 1:NB_GENOTYPES_PER_CHROM
            write(file2, "$(t);$(deme);$(chromType);$(chromosomeFrequencies[deme,chromType]) \n")
        end
    end
    TD = 500
    if (t==TD || t==TstopMigration+TD) && disassortativeMatingScript==true
        percentDIni = 0.05
        genotypeFrequencies = implementDgenotype(genotypeList, genotypeFrequencies, percentDIni)
    end
        # migration 
    if t<TstopMigration
        genotypeFrequencies = getFrequenciesAfterMigrationParapatry(genotypeList, genotypeFrequencies, genotypeOtherDemeVector, migrationRate, migrationUnidirectionnal)
    end
    
        # reproduction
    fractionMating = getMatingFraction(genotypeFrequencies, disassortativeLevel)
    genotypeFrequencies = getFrequenciesOffspringDeterministic(genotypeList, genotypeFrequencies, fractionMating, recombinationMatrix)
    
    # selection
    genotypeFrequencies = getFrequenciesAfterSelection(genotypeList, genotypeFrequencies, selectionGenotypeVector)
    
    # genetic drift
    if boolFinitePopulation==1
      genotypeFrequencies = getFrequenciesFinitePopulation(genotypeFrequencies, carryingCapacity, mutationEffectsMatrix, mutationRatePerIndividual, inversionEffectsMatrix, inversionRatePerIndividual)
    end
    
    t = t+1
  end

  close(file)
  close(file2)
  
end  


########### -----------------------------------------------


function csvSensitivityAnalysis3D(genotypeList,
                recombinationRate,
                recombinationRateRD,
                disassortativeMatingScript,
                disassortativeLevel,
                migrationRate,
                migrationUnidirectionnal,
                boolFinitePopulation,
                carryingCapacity,
                mutationRate,
                inversionRate,
                inversionType,
                percentBAIni,
                selectGoodHomozygotes,
                selectBadHomozygotes,
                selectHeterozygotes,
                additionalSelectInvHomozygotes,
                TstopMigration,
                Tend,
 	       repJulia,
	       repBash,
               nameParam1,
               nameParam2,
               nameParam3,
               nameParam4,
               nameParam5,
               valuesParam1,
               valuesParam2,
               valuesParam3,
               valuesParam4,
               valuesParam5,
               nameFile)

  if isfile("Data/SensitivityAnalysis"*nameFile*".txt")==true
    rm("Data/SensitivityAnalysis"*nameFile*".txt")
  end
  file = open("Data/SensitivityAnalysis"*nameFile*".txt", "w")
  if isfile("Data/SensitivityAnalysis2"*nameFile*".txt")==true
    rm("Data/SensitivityAnalysis2"*nameFile*".txt")
  end
  file2 = open("Data/SensitivityAnalysis2"*nameFile*".txt", "w")

  write(file,"param1;param2;param3;param4;param5;rep;state;genotype;frequency \n")
  write(file2,"param1;param2;param3;param4;param5;rep;state;deme;chromosome;frequency \n")
  
  println(nameFile)  
  
  for param1 in valuesParam1
    println(param1)
    
    if nameParam1=="selectGoodHomozygotes"
      selectGoodHomozygotes = param1
    elseif nameParam1=="selectBadHomozygotes"
      selectBadHomozygotes = param1
    elseif nameParam1=="selectHeterozygotes"
      selectHeterozygotes = param1
    elseif nameParam1=="additionalSelectInvHomozygotes"
      additionalSelectInvHomozygotes = param1
    elseif nameParam1=="recombinationRate"
      recombinationRate = param1
    elseif nameParam1=="mutationRate"
      mutationRate = param1
    elseif nameParam1=="percentBAIni"
      percentBAIni = param1
    elseif nameParam1=="migrationRate"
      migrationRate = param1
    elseif nameParam1=="inversionRate"
      inversionRate = param1
    elseif nameParam1=="disassortativeLevel"
      disassortativeLevel = param1
    elseif nameParam1=="recombinationRateRD"
      recombinationRateRD = param1
    end    
    
    for param2 in valuesParam2
        if nameParam2=="selectGoodHomozygotes"
            selectGoodHomozygotes = param2
        elseif nameParam2=="selectBadHomozygotes"
            selectBadHomozygotes = param2
        elseif nameParam2=="selectHeterozygotes"
            selectHeterozygotes = param2
        elseif nameParam2=="additionalSelectInvHomozygotes"
            additionalSelectInvHomozygotes = param2
        elseif nameParam2=="recombinationRate"
            recombinationRate = param2
        elseif nameParam2=="mutationRate"
            mutationRate = param2
        elseif nameParam2=="percentBAIni"
            percentBAIni = param2
        elseif nameParam2=="migrationRate"
            migrationRate = param2
        elseif nameParam2=="inversionRate"
            inversionRate = param2
        elseif nameParam2=="disassortativeLevel"
            disassortativeLevel = param2
        elseif nameParam2=="recombinationRateRD"
            recombinationRateRD = param2
        end
        
      for param3 in valuesParam3
        if nameParam3=="selectGoodHomozygotes"
            selectGoodHomozygotes = param3
        elseif nameParam3=="selectBadHomozygotes"
            selectBadHomozygotes = param3
        elseif nameParam3=="selectHeterozygotes"
            selectHeterozygotes = param3
        elseif nameParam3=="additionalSelectInvHomozygotes"
            additionalSelectInvHomozygotes = param3
        elseif nameParam3=="recombinationRate"
            recombinationRate = param3
        elseif nameParam3=="mutationRate"
            mutationRate = param3
        elseif nameParam3=="percentBAIni"
            percentBAIni = param3
        elseif nameParam3=="migrationRate"
            migrationRate = param3
        elseif nameParam3=="inversionRate"
            inversionRate = param3
        elseif nameParam3=="disassortativeLevel"
            disassortativeLevel = param3
        elseif nameParam3=="recombinationRateRD"
            recombinationRateRD = param3
        end    
	      for param4 in valuesParam4
	        if nameParam4=="selectGoodHomozygotes"
	            selectGoodHomozygotes = param4
	        elseif nameParam4=="selectBadHomozygotes"
	            selectBadHomozygotes = param4
	        elseif nameParam4=="selectHeterozygotes"
	            selectHeterozygotes = param4
	        elseif nameParam4=="additionalSelectInvHomozygotes"
	            additionalSelectInvHomozygotes = param4
	        elseif nameParam4=="recombinationRate"
	            recombinationRate = param4
	        elseif nameParam4=="mutationRate"
	            mutationRate = param4
	        elseif nameParam4=="percentBAIni"
	            percentBAIni = param4
	        elseif nameParam4=="migrationRate"
	            migrationRate = param4
	        elseif nameParam4=="inversionRate"
	            inversionRate = param4
 		elseif nameParam4=="disassortativeLevel"
	            disassortativeLevel = param4
	        elseif nameParam4=="recombinationRateRD"
	            recombinationRateRD = param4
	        end 

			for param5 in valuesParam5
		        if nameParam5=="selectGoodHomozygotes"
		            selectGoodHomozygotes = param5
		        elseif nameParam5=="selectBadHomozygotes"
		            selectBadHomozygotes = param5
		        elseif nameParam5=="selectHeterozygotes"
		            selectHeterozygotes = param5
		        elseif nameParam5=="additionalSelectInvHomozygotes"
		            additionalSelectInvHomozygotes = param5
		        elseif nameParam5=="recombinationRate"
		            recombinationRate = param5
		        elseif nameParam5=="mutationRate"
		            mutationRate = param5
		        elseif nameParam5=="percentBAIni"
		            percentBAIni = param5
		        elseif nameParam5=="migrationRate"
		            migrationRate = param5
		        elseif nameParam5=="inversionRate"
		            inversionRate = param5
		        elseif nameParam5=="disassortativeLevel"
		            disassortativeLevel = param5
		        elseif nameParam5=="recombinationRateRD"
		            recombinationRateRD = param5
		        end 
  
			for rep in 1:repJulia
			      repAll = (repBash-1)*repJulia+rep

			      #mutation effects
			      mutationEffectsMatrix = createMutationEffectsMatrix(genotypeList, mutationRate)
			      mutationRatePerIndividual = calculateMutationRatePerIndividual(genotypeList, mutationRate)   

			      #inversion effects
			      if inversionType==0
				inversionEffectsMatrix = createInversionEffectsMatrixTwoWays(genotypeList, inversionRate)
				inversionRatePerIndividual = calculateInversionRatePerIndividualTwoWays(genotypeList, inversionRate) 
			      elseif inversionType==1
				inversionEffectsMatrix = createInversionEffectsMatrixOneWay(genotypeList, inversionRate)
				inversionRatePerIndividual = calculateInversionRatePerIndividualOneWay(genotypeList, inversionRate) 
			      end
			      #recombination rate
			      recombinationMatrix = createRecombinationMatrix(genotypeList, recombinationRate, recombinationRateRD)

			      #selection
			      selectionGenotypeVector = createSelectionGenotypeVector(selectGoodHomozygotes, selectBadHomozygotes, selectHeterozygotes, additionalSelectInvHomozygotes)

			      #initialization
			      genotypeFrequencies = initializeGenotypeFrequencies(genotypeList, typeOfInitialization, percentBAIni, recombinationMatrix, selectionGenotypeVector, migrationRate, migrationUnidirectionnal)

			      # with migration
			      genotypeFrequencies = freqAfterSimulation(genotypeList, genotypeFrequencies, recombinationMatrix, mutationEffectsMatrix, mutationRatePerIndividual, inversionEffectsMatrix, inversionRatePerIndividual,  selectionGenotypeVector, migrationRate, migrationUnidirectionnal, TstopMigration)

			      for genotype in 1:NB_GENOTYPES
				write(file, "$(param1);$(param2);$(param3);$(param4);$(param5);$(repAll);0;$(genotype);$(genotypeFrequencies[genotype]) \n")
			      end
			      chromosomeFrequencies = getChromosomeFrequencyEachDeme(genotypeFrequencies, genotypeList)
			      for deme in 1:NB_DEMES
				  for chromType in 1:NB_GENOTYPES_PER_CHROM
				    write(file2, "$(param1);$(param2);$(param3);$(param4);$(param5);$(repAll);0;$(deme);$(chromType);$(chromosomeFrequencies[deme,chromType]) \n")
				  end
			      end      
			      
			      # without migration
			      newMigrationRate = 0.0
			      genotypeFrequencies = freqAfterSimulation(genotypeList, genotypeFrequencies, recombinationMatrix, mutationEffectsMatrix, mutationRatePerIndividual, inversionEffectsMatrix, inversionRatePerIndividual, selectionGenotypeVector, newMigrationRate, migrationUnidirectionnal, Tend-TstopMigration)
			      
			      for genotype in 1:NB_GENOTYPES
				write(file, "$(param1);$(param2);$(param3);$(param4);$(param5);$(repAll);1;$(genotype);$(genotypeFrequencies[genotype]) \n")
			      end
			      chromosomeFrequencies = getChromosomeFrequencyEachDeme(genotypeFrequencies, genotypeList)
			      for deme in 1:NB_DEMES
				  for chromType in 1:NB_GENOTYPES_PER_CHROM
				    write(file2, "$(param1);$(param2);$(param3);$(param4);$(param5);$(repAll);1;$(deme);$(chromType);$(chromosomeFrequencies[deme,chromType]) \n")
				  end
			      end     
			end
		    end
		end
	      end
	  end  
  end

  close(file)
  close(file2)
  
end  


########### ---------------------------------------------------------------------
