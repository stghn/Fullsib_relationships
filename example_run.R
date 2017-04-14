grm='../AGDFieldGenomics/GMAT/MH2015GCMS/MH2015CMS_asreml.grm'
genolinkped='../AGDFieldGenomics/GMAT/MH2015GCMS/grm.orig_codedIDs'
pedigree='../AGDFieldGenomics/GMAT/MH2015GCMS/pedigreeMH2015G.txt'

outname='AGDMH2015G'

source('extfullsib_GRM.R')
AGDMH <- extfullsib_GRM(grm ,pedigree,genolinkped,exclfamsize=13,outname)
