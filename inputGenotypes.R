NbGenotypes = 128 #ModifPJ
NbGenotypesPerChrom = 8#ModifPJ
NbDemes = 2#ModifPJ

### list of genotypes: pair of chromosomes
genotypeList =  array(0,dim = c(NbGenotypes,3))
genotypeNb = 1
for(deme in 1:NbDemes){
  for(i in 1:NbGenotypesPerChrom){
    for(j in 1:NbGenotypesPerChrom){
      genotypeList[genotypeNb,1] = i 
      genotypeList[genotypeNb,2] = j 
      genotypeList[genotypeNb,3] = deme 
      genotypeNb = genotypeNb + 1
    } 
  }
}

deme2list = c();
chromBAHomlist = c();
chromBAHetlist = c();
chromABHomlist = c();
chromABHetlist = c();
chromabHomlist = c();
chromOtherlist = c()
chrom1=1;chrom2=2;demeIndex=3;
alleleBA=1;alleleAB=5;alleleab=8;alleleba=4#ModifPJ
for(genotype in 1:NbGenotypes){
  if(genotypeList[genotype,demeIndex]==2){
    deme2list = c(deme2list,genotype);
  }
  if(genotypeList[genotype,chrom1]==alleleBA & genotypeList[genotype,chrom2]==alleleBA){ chromBAHomlist = c(chromBAHomlist,genotype) }
else if(genotypeList[genotype,chrom1]==alleleAB & genotypeList[genotype,chrom2]==alleleAB){ chromABHomlist = c(chromABHomlist,genotype) }
else if(genotypeList[genotype,chrom1]==alleleab & genotypeList[genotype,chrom2]==alleleab){ chromabHomlist = c(chromabHomlist,genotype) }
else if(genotypeList[genotype,chrom1]==alleleBA | genotypeList[genotype,chrom2]==alleleBA){ chromBAHetlist = c(chromBAHetlist,genotype) }
  # }else if(genotypeList[genotype,chrom1]==alleleAB | 
  #          genotypeList[genotype,chrom2]==alleleAB)
  #   chromABHetlist = c(chromABHetlist,genotype)    
  
else if(genotypeList[genotype,chrom1]==alleleAB | 
           genotypeList[genotype,chrom2]==alleleAB){
    chromABHetlist = c(chromABHetlist,genotype)
  }else{
    chromOtherlist = c(chromOtherlist,genotype)
  }
}


### list of chromosomes for SingleTimeSeries2

chromBAlist2 = c(1);#ModifPJ
chrombalist2 = c(4);#ModifPJ
chromInvOthlist2 = c(2,3);#ModifPJ
chromABlist2 = c(5);#ModifPJ
chromablist2 = c(8);#ModifPJ
chromOtherlist2 = c(6,7)#ModifPJ
