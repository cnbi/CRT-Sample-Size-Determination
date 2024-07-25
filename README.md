# Power analysis in cluster randomized trials
This is a research compendium containing all the scripts, manuscript, and illustrations used in the project “Sample Size Determination for Cluster Randomised Trials with the Bayes Factor”. The aim of this project was to develop a method and tools for sample size determination of cluster randomised trials when using the Approximated Adjusted Fractional Bayes Factor for hypothesis testing. 

## Repository Structure
The following graph displays the organization of the files and folders related to this research study:
-	**docs**
      -	**manuscript** – manuscript source and biblatex with all the references
-	**results**
      -	**figures** – Illustrations used in the manuscript.
      -	**simulation results** – Data sets with simulation results.
-	**scripts**
      -	Functions
      -	Scripts with simulations
      -	**plots** – Scripts used to create plots for presentations.
      -	**performance** – Scripts for profiling the functions created.
-	CITATION.md - 
-	LICENSE
-	README.md
-	CITATION.cff - Citation metadata for this repository citation.

## R package
This project is part of a bigger project named “Bayesian sample size calculation for trials with multilevel data”. All the functions developed to determine the sample size for cluster randomised trials are available in the R package [SSD_Bayes_ML]( https://github.com/ulrichlosener/SSD_Bayes_ML)

## Shiny App
In the manuscript can be seen a fraction of the simulation results given the limitations in space. For this reason, we have created a Shiny app in which is possible to explore the effects of different combinations of the intraclass correlation, Bayes factor threshold, cluster size, number of clusters, effect sizes and fractions of information (b).
[Bayesian Sample Size Determination: Cluster Randomised Trials](https://utrecht-university.shinyapps.io/BayesSamplSizeDet-CRT/)

## Citing
Citing repository:

Citing article:
**BibTex**
@misc{barragan-ibanez_moerbeek_2024,
 title={Method for Sample Size Determination for Cluster-Randomized Trials Using the Bayes Factor},
 url={osf.io/preprints/psyarxiv/kcz7v},
 DOI={10.31234/osf.io/kcz7v},
 publisher={PsyArXiv},
 author={Barragan-Ibanez, Camila N and Moerbeek, Mirjam},
 year={2024},
 month={Jul}
}

**APA**
Barragan-Ibanez, C. N., & Moerbeek, M. (2024, July 7). Method for Sample Size Determination for Cluster-Randomized Trials Using the Bayes Factor. https://doi.org/10.31234/osf.io/kcz7v

## License
This project is licensed under the terms of the [GNU GPLv3](/LICENSE)
