# Supplementary Material for ["Dynamics of a plant-pollinator network: Extending the Bianconi-Barabási Model"](https://rdcu.be/dMeuQ) 

This repository contains the supplementary material to the paper "Dynamics of a plant-pollinator network: Extending the Bianconi-Barabási Model" by William J. Castillo, Laura A. Burkle and Carsten F. Dormann, published in the Special Issue of the 12th International Conference on Complex Networks and their Applications of the journal Applied Network Science, DOI: 10.1007/s41109-024-00636-0.

## Contents

- **Data**: Datasets used in the analysis.
- **Code**: Scripts and functions for data processing.

## Repository Structure


### Data

The information is ditributed in the next .CSV files: 
  - `bee_interactions` A table containing the interactions between bee and plant individuals during the summer of 2014.
  - `BeeTraits` Contains data on each bee species, detailing the average intertegular distances, average bee size, and species abundance.
  - `PlantTraits` Includes data on each plant species, detailing the average height (from ground level to canopy), the average number of flowers per stem, and the average width of the flower head.

Additionally, `ListMatrixInteractions.rdata` is an R list where each element is a matrix representing the number of interactions at each time step. The provided functions process the network information using these matrices, which are based on the data from `bee_interactions`. 

### Code

The script `funktionen.R` contains the functions used in the data analysis, including those for generating the degree and strength dynamics of the network assembly.

If you have any comments or questions, please don't hesitate to contact us at: <william.castillo@biom.uni-freiburg.de>
