# Supplementary Material for "Dynamics of a plant-pollinator network: Extending the Bianconi-Barabási Model"

This repository contains the supplementary material for the paper "Dynamics of a plant-pollinator network: Extending the Bianconi-Barabási Model" by William J Castillo, Laura A. Burkle, and Carsten F Dormann, submitted to the "Special Issue of the 12th International Conference on Complex Networks and their Applications".

## Contents

- **Data**: Datasets used in the analysis.
- **Code**: Scripts and functions for data processing.

## Repository Structure


### Data

The information is ditributed in the next .CSV files: 
  - `bee_interactions` A table containing interactions between bee and plant species during the summer of 2014.
  - `BeeTraits` Contains the intertegular distances, bee size, and abundance of each bee species.
  - `PlantTraits` Contains the height from the ground level to the canopy of plants, the number of flowers per stem, and the width of the flower head.

Each group of functional traits is divided by the number of species.

### Code

The script contains functions written in R to generate the time series representing the network assembly. 
