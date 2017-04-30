grm='../AGDFieldGenomics/GMAT/MHgeno/AGD_GS_asreml.grm'
genolinkped='../AGDFieldGenomics/GMAT/MHgeno/grm.orig_codedIDs'
pedigree='../AGDFieldGenomics/GMAT/MHgeno/pedigree.txt'
outname='G2015_HS'

source('extractFHsibrel_GRM.R')
AGDMH <- exFHsibrel_GRM(grm,pedigree,genolinkped,typefam = 'HS',exclfamsize=1,outname)


grm='../AGDFieldGenomics/GMAT/MHgeno/AGD_GS_asreml.grm'
genolinkped='../AGDFieldGenomics/GMAT/MHgeno/grm.orig_codedIDs'
pedigree='../AGDFieldGenomics/GMAT/MHgeno/pedigree.txt'
outname='G2015_FS'

source('extractFHsibrel_GRM.R')
AGDMH <- exFHsibrel_GRM(grm,pedigree,genolinkped,typefam='FS',exclfamsize=1,outname)
