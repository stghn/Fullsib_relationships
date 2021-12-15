## Extracting Fullsib and Halfsib relationships from a G matrix and plot the relationships 

### The script requires the following
        
  - **grm**: Genomic relationships in long column and rowwise **(triplet)** format -- **not in matrix format** 
  - **pedigree**: Pedigree file **(ID, sire and Dam, missing parent = NA)**
  - **genolinkped**: Linking animals in the G-matrix to animals in pedigree **(2 columns IID-pedigree, IID-G matrix)**
  - **typefam**: Type of familial relationships to extract **(only two types are supported *fullsib ='FS'* and *halfsib='HS'*)**
  - **exclfamsize**: Exclude families with **n number of animals** (e.g. a family with 1 animal)
  - **outname**: Output name of the graph and the output file      
        
### Output
  - A graph with the relationships
  - A file contaning the extract full/half sib relationships
---
### Example run

```R
source('extractFHsibrel_GRM.R')  ## Executing source function to load R script
example_HS <- exFHsibrel_GRM(ggrm='example.grm',genolinkped='orig_recodedIDs.txt',
                            pedigree='pedigree.txt',typefam='FS',exclfamsize=1,outname='G2015_FS')
```
 
