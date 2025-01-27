# salvia_develop1
This repository includes the raw data from "Ontogenetic mechanisms of differentiation in two Salvia species with different pollinators" by A. Davies &amp; S. Benitez-Vieyra, alongwith the R script to replicate all the analyses.

File and variable names.   
   
* celulas_num.csv. This file contains information about the relation between flower size and cell number in two *Salvia* species.   
   + species. Plant species: S.guaranitica, *Salvia guaranitica*; S.stachydifolia, *S. stachydifolia*   
   + ID. Individual ID number.
   + bud. Bud or flower (nested in ID).
   + CN. Cell number.
   + TL. total length (from the base to the insertion point of the stamens).   
   
* celulas_tam.csv. This file contains information about the relation between flower size and cell size in two *Salvia* species.    
   + species. Plant species: S.guaranitica, *Salvia guaranitica*; S.stachydifolia, *S. stachydifolia*  
   + ID. Individual ID number.
   + bud. Bud or flower (nested in ID).
   + cell. Cell ID (nested in bud).
   + region. Either basal, medium or distal part of the flower tube.
   + W. Cell width.
   + L. Cell length.
   + TL. total length (from the base to the insertion point of the stamens).  

* trayetorias_ont.csv. This file contains information about the relation between corolla tube and corolla upper lip size in two *Salvia* species.     
   + species. Plant species: gua, *Salvia guaranitica*; sta, *S. stachydifolia*   
   + ID. Individual ID number.
   + bud. Bud or flower (nested in ID).
   + CTL. Corolla tube length.
   + ULL. Corolla upper lip length.   
   + TL. total length

* calyx_corolla.csv. This file contains information about the relation between corolla tube length, corolla length and calyx length in 30 flower buds of two *Salvia* species.
   + id. Individual ID number.   
   + sp. S.guaranitica, *Salvia guaranitica*; S.stachydifolia, *S. stachydifolia*  
   + calyx_length. In mm.
   + corolla_length. In mm.   
   + tube_length. In mm.
 
* crecimiento.csv. This file contains information about the relation between flower size and days to anthesis in two *Salvia* species.
   + sp. S.guaranitica, *Salvia guaranitica*; S.stachydifolia, *S. stachydifolia*
   + ID. Individual bud.
   + plant. Individual plant.
   + day_orig. Day of measurement in incremental sequence.
   + long_calix. Calyx length in mm.
   + long_tube. Corolla tube length in mm. Only for flowers.   
   + long_tot. Total corolla length in mm.
   + cat. Factor. b = bud, f = flower.
   + day_pre. Days to anthesis.
   + hours. Hours to anthesis.

* script_rev.R. This file contains all the details needed to replicate statistical analyses and graphs.   
