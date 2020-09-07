# Scripts for various population genomic analyses

### Population genetics summary statistics
  * F<sub>ST</sub> - `popgenome_statistics_calculations.R`
  * D<sub>XY</sub> - `popgenome_statistics_calculations.R`
  * Nucleotide diversity (π) - `popgenome_statistics_calculations.R`
  * Fis - `calculate_Fis.R`

D<sub>a</sub> was calculated by subtracting D<sub>XY</sub> from female pi.

### Liftover to curated chr12 positions
  * A liftover function (`update_chr12_liftover.R`) is provided to update chromosome 12 positions from published assembly to curated positions. 
