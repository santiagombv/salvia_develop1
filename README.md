# salvia_develop1
This repository includes the raw data from "Ontogenetic mechanisms of differentiation in two Salvia species with different pollinators" by A. Davies &amp; S. Benitez-Vieyra, alongwith the R script to replicate all the analyses.

File and variable names.   
   
* celulas_num.csv. This file contain information about the relation between flower size and cell number in two *Salvia* species.   
   + species. Plant species: guar, *Salvia guaranitica*; stach, *S. stachydifolia*   
   + ID. Individual ID number.
   + bud. Bud or flower (nested in ID).
   + CN. Cell number.
   + TL. total length (from the base to the insertion point of the stamens).   
   
* celulas_tam.csv. This file contain information about the relation between flower size and cell size in two *Salvia* species.    
   + species. Plant species: *S. guaranitica*, *S. stachydifolia*   
   + ID. Individual ID number.
   + bud. Bud or flower (nested in ID).
   + cell. Cell ID (nested in bud).
   + region. Either basal, medium or dital part of the flower tube.
   + W. Cell width.
   + L. Cell length.
   + TL. total length (from the base to the insertion point of the stamens).  

* trayetorias_ont.csv. This file contain information about the relation between corolla tube and corolla upper lip size in two *Salvia* species.     
   + species. Plant species: gua, *Salvia guaranitica*; sta, *S. stachydifolia*   
   + ID. Individual ID number.
   + bud. Bud or flower (nested in ID).
   + CTL. Corolla tube length.
   + ULL. Corolla upper lip length.   
   + TL. total length

* script.R. This file contains all the details needed to replicate statistical analyses and graphs.   
