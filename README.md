![Striking-image](https://github.com/mirzaei-s/microbial_interactions/assets/150903671/98da836b-e0b4-41c7-b7af-e7f2f984c45b)

# microbial_piarwise_interactions
This function explores which commonly consumed metabolites two species have the potential to compete for. the results represnt as competition rate.
Also, check the possibility of parasitic and commensal interactions for each given-consumer metabolites between two species.

% INPUT:

 model:       A community COBRA model structure with the following fields

                    S - Stoichiometric matrix
                    b - Right hand side
                    c - Objective coefficients
                    lb - Lower bounds
                    ub - Upper bounds
                    
 % OUTPUTS:
 
         CompeteList: List of metabolites  which two species compete.
         Commensal_species1 : List of metabolites that species2 demans from species1. ("1 can provide (commensal), 0 cannot ")
         Commensal_species2 : List of metabolites that species1 demans from species2. 

                    

                    
[CompeteList, Commensal_species1,Commensal_species2]  = PairwiseInteraction(model1,model2); Call function with two model.

# Citation
This repository contains the code and configuration files reproducing the study cases described in (please, cite us if you use it in your work):

    Mirzaei, Soraya, and Mojtaba Tefagh. "GEM-based computational modeling for exploring metabolic interactions in a microbial community." PLOS Computational Biology 20.6 (2024): e1012233.

