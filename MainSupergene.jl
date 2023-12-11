include("./AnalysisModelSupergene.jl")



########### PARAMETERS --------------------------------------------------------------------- #Modify this for changing the parameter of each simulation


const typeOfInitialization = "OnlyABab" # Uniform, UniformNoBA, OnlyAB, Onlyab, OnlyBA, OnlyABab or InvasionBA
const percentBAIni = 0.05  # useless

const disassortativeMatingScript = false
const migrationUnidirectionnal = false

const boolFinitePopulation = 1
carryingCapacity = 10000

mutationRate = 1e-6
inversionRate = 1e-6
const inversionType = 1  # 0: both ways; 1:one way

const disassortativeLevel = 0.0
const recombinationRateRD = 0.5

migrationRate = 0.05
selectGoodHomozygotes = 1.0
selectBadHomozygotes = 0.5
selectHeterozygotes = 0.2
additionalSelectInvHomozygotes = 0.95
recombinationRate = 0.2

########### DO NOT TOUCH  ---------------------------------------------------------------------


if disassortativeMatingScript==true
    include("./FunctionsSupergeneD.jl")
else
    include("./FunctionsSupergene.jl")
end
const genotypeList = createGenotypes()
const phenotypeList = createPhenotypesEco()
const phenotypeFardList = createPhenotypesFardeau()
recombinationMatrix = createRecombinationMatrix(genotypeList, recombinationRate, recombinationRateRD)
selectionGenotypeVector = createSelectionGenotypeVector(selectGoodHomozygotes, selectBadHomozygotes, selectHeterozygotes, additionalSelectInvHomozygotes)
const genotypeOtherDemeVector = createVectorGenotypeOtherDeme(genotypeList)



########### SIMULATIONS ---------------------------------------------------------------------


### ONE SIMULATION ### Uncomment this to perform a single simulation

#  const nameFile = "test7one"
#  const TstopMigration = 5000
#   const Tend = 10000
#  @time csvOneSimulation(genotypeList, recombinationMatrix, disassortativeMatingScript, disassortativeLevel, selectionGenotypeVector, genotypeOtherDemeVector, migrationRate, migrationUnidirectionnal, boolFinitePopulation, carryingCapacity, mutationRate, inversionRate, inversionType, TstopMigration, Tend, nameFile)
#
### SENSITIVITY ANALYSIS #Perform analyses exploring a parameter space. Modify the parameter arrays to change the parameter space explored

const nameFile = ARGS[1]

const TstopMigration = 5000
const Tend = 10000
const repJulia = parse(Int32, ARGS[2])
const repBash = parse(Int32, ARGS[3])

const nameParam1 = "migrationRate"
const valuesParam1 = [0.05 0.1 0.15] 

const nameParam2 = "recombinationRate"
const valuesParam2 = [0.0 0.1 0.2 0.3 0.4 0.5]

const nameParam3 = "additionalSelectInvHomozygotes" #Mutation load
const valuesParam3 = [0.0 0.05 0.1 0.15 0.20 0.25 0.3 0.35]

const nameParam4 = "selectHeterozygotes" #Fitness of recombinant
const valuesParam4 = [0.2]

const nameParam5 = "selectBadHomozygotes" #Fitness of immigrant genotypes (e.g. ab/ab in pop 1, AB/AB in pop 2)
const valuesParam5 = [0.5]

csvSensitivityAnalysis3D(genotypeList,
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
