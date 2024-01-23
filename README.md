# microbial_piarwise_interactions
This function explores which commonly consumed metabolites two species have the potential to compete for. the results represnt as competition rate.
Also, check the possibility of parasitic and commensal interactions for each given-consumer metabolites between two species.

% INPUT:

%    model:       A community COBRA model structure with the following fields

                    S - Stoichiometric matrix
                    b - Right hand side
                    c - Objective coefficients
                    
                    lb - Lower bounds
                    
                    ub - Upper bounds
                    
Call PairwiseInteraction() function with two model.
