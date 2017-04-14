## Extracting Fullsibs from a G matrix and plotting the results

### The script requires the following
        
        - grm             ==> Genomic relationships in long column and rowwise format (not in matrix format) 
        - pedigree        ==> pedigree file (ID, sire and Dam, missing parent = NA)
        - genolinkped     ==> linking animals in the G matrix to pedigre  (2 columns IID-pedigree, IID-G matrix)
        - exclfamsize     ==> exclude families with n number of animals (e.g. a family with 1 animal)
        - outname         ==> output name of the graph and the output file 
        
        
### Output
        - A graph with the relationships
        - A file contaning the extract full/half sib relationships
