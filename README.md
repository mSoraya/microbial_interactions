![Striking-image](https://github.com/mirzaei-s/microbial_interactions/assets/150903671/fd1ae370-aebc-432e-a7ab-1b01d07fc4b3)
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
