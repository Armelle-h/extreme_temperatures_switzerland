# Extreme spatial temperature events over Switzerland

---

## About the Project

### Description

This repository contains the code associated with my Master thesis : "Extreme spatial temperature events over Switzerland", completed as a part of my academic studies as an Applied Mathematics Master student at EPFL.

Predicting high temperatures is crucial for risk management, as they can impact public health, agriculture, infrastructure, and so on by causing health risks such as heatstroke, stress on crops and livestock, and damage to roads and bridges. The thesis consists of investigate the evolution of extreme temperatures in Switzerland between 1971 and 2022 by adapting the methodology described in the paper the "Inference for extreme spatial temperature events in a changing climate with application to Ireland" by Healy et al. 


(to remove: We model the tail of the marginal distribution with the generalized Pareto distribution and the extreme temperatures dependencies with a $r$-Pareto process.)

---

## Code:
- The folder [marginal_model] contains ...
- The folder [spatial_model] contains ..

## Usage 
- The folder [marginal_model] should be run before the folder [spatial_model].
- In the folder [marginal_model], the files related to the bulk model should be run before the files related to the gpd model.
- In the folder [spatial_model], (to complete once I've done file cleaning)
- The file [simulPareto.R] in the folder [mvPot] is an adaptation of the function "simulPareto" from the package MvPot (https://github.com/r-fndv/mvPot)

## Dependencies:
The main packages used are 
- tidyverse
- evgam
- (to complete)
  
---
## Main references
- reference 1 

## Acknowledgments

- My sincere gratitude to Professor Anthony C. Davison, my master thesis advisor and supervisor at EPFL.
- I am grateful to Dr. Daire Healy for providing public access to their code (https://github.com/dairer/Extreme-Irish-Summer-Temperatures/tree/main), associated with their paper "Inference for extreme spatial temperature events in a changing climate with application to Ireland" which highly supported this work.
