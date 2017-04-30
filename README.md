## Extracting Fullsib and Halfsib relationships from a G matrix and plotting the results


### The script requires the following
        
        - grm             ==> Genomic relationships in long column and rowwise format (not in matrix format) 
        - pedigree        ==> pedigree file (ID, sire and Dam, missing parent = NA)
        - genolinkped     ==> linking animals in the G matrix to pedigre  (2 columns IID-pedigree, IID-G matrix)
        - typefam         ==> type of familial relationships to extract (only two types are supported fullsib =FS' and hlfsib='HS' )
        - exclfamsize     ==> exclude families with n number of animals (e.g. a family with 1 animal)
        - outname         ==> output name of the graph and the output file 
        
        
### Output
        - A graph with the relationships
        - A file contaning the extract full/half sib relationships


### example run

      source the script 
        source('extractFHsibrel_GRM.R')
        
        example_HS <- exFHsibrel_GRM(ggrm='example.grm',genolinkped='orig_recodedIDs.txt',
        pedigree='pedigree.txt',typefam='FS',exclfamsize=1,outname='G2015_FS')
