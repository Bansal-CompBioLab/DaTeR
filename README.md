# DaTeR: Phylogenetic dating using relative time constraints

## **Description**
DaTeR (short for “Dating Trees using Relative constraints”) is a program for improved dating of microbial species phylogenies using *relative* time constraints (e.g., obtained from high-confidence horizontal gene transfer events). Traditional phylogenetic dating approaches make use of *absolute* time constraints, which provide lower and/or upper bounds for one or more nodes of the underlying phylogeny, but are unable to use relative constraints that specify that some node *x* must be at dated to be at least as old as some other node *y*. DaTeR takes as input a collection of chronograms sampled from the posterior using any traditional Bayesian phylogenetic dating approach (based on only absolute time calibrations), along with a set of curated relative time constraints, and minimally error-corrects each input chronogram to ensure compatibility with all available relative time constraints. It then outputs the individual error-corrected chronogram samples as well as an aggregated, final chronogram. DaTeR uses a constrained optimization framework and computes a minimal deviation from assigned node dates or branch lengths (representing time) under three appropriately designed candidate objective functions. Further technical details appear in the paper cited below. DaTeR was implemented by Abhijit Mondal and is available open-source under GNU GPL.

This repository includes complete source code, user manual, and test data. For reference, a copy of the user manual is available from this GitHub link: [DaTeR user manual]() 


## **Citation information**
DaTeR can be cited as follows:

<a href="https://doi.org/10.1093/bioinformatics/btad084">DaTeR: Error-Correcting Phylogenetic Chronograms Using Relative Time Constraints</a><br>
Abhijit Mondal, L. Thiberio Rangel, Jack Payette, Gregory P. Fournier, Mukul S. Bansal<br>
*Bioinformatics*, 39(2), btad084, 2023
