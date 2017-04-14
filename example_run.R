grm='../AGDFieldGenomics/GMAT/MHgeno/AGD_GS_asreml.grm'
genolinkped='../AGDFieldGenomics/GMAT/MHgeno/grm.orig_codedIDs'
pedigree='../AGDFieldGenomics/GMAT/MHgeno/pedigree.txt'
outname='G2015_FS'

source('extfullsib_GRM.R')
AGDMH <- extfullsib_GRM(grm ,pedigree,genolinkped,exclfamsize=1,outname)
