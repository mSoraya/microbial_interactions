function []  = PairwiseInteraction(model1,model2)

% INPUT:
%    model:            COBRA model structure
%

% OUTPUTS:
[model1_Irrev, ~, ~, ~] = convertToIrreversible(model1);
[model2_Irrev, ~, ~, ~] = convertToIrreversible(model2);


[BlockedRxns] = identifyFastBlockedRxns(model1_Irrev,model1_Irrev.rxns);
index=find(ismember(model1_Irrev.rxns,BlockedRxns));
model1_Irrev.S(:,index)=[];
model1_Irrev.rxns(index)=[];
model1_Irrev.c(index)=[];
model1_Irrev.lb(index)=[];
model1_Irrev.ub(index)=[];

clear BlockedRxns index;
[BlockedRxns] = identifyFastBlockedRxns(model2_Irrev,model2_Irrev.rxns);
index=find(ismember(model2_Irrev.rxns,BlockedRxns));
model2_Irrev.S(:,index)=[];
model2_Irrev.rxns(index)=[];
model2_Irrev.c(index)=[];
model2_Irrev.lb(index)=[];
model2_Irrev.ub(index)=[];


%%  Preprocesing
[minFlux, maxFlux] = fluxVariability(model1_Irrev);
model1_Irrev.minFlux=minFlux;
model1_Irrev.maxFlux=maxFlux;

clear minFlux maxFlux;
[minFlux, maxFlux] = fluxVariability(model2_Irrev);
model2_Irrev.minFlux=minFlux;
model2_Irrev.maxFlux=maxFlux;

biomass=find(model1_Irrev.c~=0);
model1_Irrev.rxns(biomass)=cellstr("Biomass");

biomass=find(model2_Irrev.c~=0);
model2_Irrev.rxns(biomass)=cellstr("Biomass");

clear minFlux maxFlux biomass;
%% Merge Two Models

model=Merge_2_Models(model1_Irrev,model2_Irrev);
[m,n]=size(model.S);
%% Competition
FBAsolution1 = optimizeCbModel(model1_Irrev,'max');
biomass1=find(contains(model.rxns,'Biomass_species1'));
model.lb(biomass1,1)=0.7*FBAsolution1.f;

FBAsolution2 = optimizeCbModel(model2_Irrev,'max');
biomass2=find(contains(model.rxns,'Biomass_species2'));
model.lb(biomass2,1)=0.7*FBAsolution2.f;

k=1;
for i=1:length(model.shared_rxn(:,1))
    temp=find(model.shared_rxn(i,:)>0);
    shared_rxn=model.shared_rxn(i,temp);
  
    Compete_rxns=model.shared_rxn(i,find(model.S( model.shared_met(i,1),shared_rxn)<0));

    if(length(Compete_rxns)==2)
        Iscompete=Competition(Compete_rxns(1,1),Compete_rxns(1,2),model);
        Competition_result{k,1}=model.rxns(Compete_rxns(1,1));
        Competition_result{k,2}=model.rxns(Compete_rxns(1,2));
        Competition_result{k,3}=Iscompete;
        Competition_result{k,4}=model.minFlux(Compete_rxns(1,1));
        Competition_result{k,5}=model.minFlux(Compete_rxns(1,2));
        Competition_result{k,6}=model.maxFlux(Compete_rxns(1,1));
        Competition_result{k,7}=model.maxFlux(Compete_rxns(1,2));
        k=k+1;
       
    end
end

iscompete = cellfun( @(X) X ==1 , Competition_result(:,3), 'UniformOutput', true );
numCompete = sum(iscompete);
[all,~]=size(Competition_result);
rate_of_Competion=numCompete /length(model.shared_rxn(:,1));
 fprintf('The rate of competiton is: %d\n',rate_of_Competion);
%% Parasitisim 


k1=1;k2=1;;all1=0;all2=0;
% for Parasitisim of species1
for i=1:length(model.shared_rxn(:,1))
    %shared_rxn=model.shared_rxn(i,:);
    temp=find(model.shared_rxn(i,:)>0);
    shared_rxn=model.shared_rxn(i,temp);
  
    rxn_species1=find(contains(model.rxns(shared_rxn),'_species1'));
    rxn_species2=find(contains(model.rxns(shared_rxn),'_species2'));

    if(any(model.S(model.shared_met(i,1),shared_rxn(1,rxn_species1)) < 0))
        if(any(model.S(model.shared_met(i,1),shared_rxn(1,rxn_species2)) > 0))
            
            s1_need=shared_rxn(1,rxn_species1((find(model.S(model.shared_met(i,1),shared_rxn(1,rxn_species1)) < 0))));
            all1=all1+1;
             s2_force=shared_rxn(1,rxn_species2((find(model.S(model.shared_met(i,1),shared_rxn(1,rxn_species2)) > 0))));
                 
            if(model.minFlux(s1_need)~=0 & model.minFlux(s1_need) > model.minFlux(s2_force))
                
                 Isparasit= Parasitisim(s2_force,biomass2,model.minFlux(s1_need),FBAsolution2.f,model);
               Species1_parasit{k1,1}=model.mets(model.shared_met(i,1));
               Species1_parasit{k1,2}=Isparasit;
               Species1_parasit{k1,3}=model.minFlux(s1_need);
               Species1_parasit{k1,4}=model.minFlux(s2_force); 
               Species1_parasit{k1,5}=model.rxns(s2_force);
               Species1_parasit{k1,6}=model.rxns(s1_need);
               
               k1=k1+1;
            end
        end
    end


    % for Parasitisim of species2
    if(any(model.S(model.shared_met(i,1),shared_rxn(1,rxn_species2)) < 0))
        if(any(model.S(model.shared_met(i,1),shared_rxn(1,rxn_species1)) > 0))
            
            s2_need=shared_rxn(1,rxn_species2((find(model.S(model.shared_met(i,1),shared_rxn(1,rxn_species2)) < 0))));
         all2=all2+1;
          s1_force=shared_rxn(1,rxn_species1((find(model.S(model.shared_met(i,1),shared_rxn(1,rxn_species1)) > 0))));
            if(model.minFlux(s2_need)~=0 & model.minFlux(s2_need) > model.minFlux(s1_force))
                model.rxns(shared_rxn(1,rxn_species2));
                 model.rxns(shared_rxn(1,rxn_species1));
               Isparasit= Parasitisim(s1_force,biomass1,model.minFlux(s2_need),FBAsolution1.f,model);
               Species2_parasit{k2,1}=model.mets(model.shared_met(i,1));
               Species2_parasit{k2,2}=Isparasit;
               Species2_parasit{k2,3}=model.minFlux(s2_need);
               Species2_parasit{k2,4}=model.minFlux(s1_force);
               Species2_parasit{k2,5}=model.rxns(s1_force);
               Species2_parasit{k2,6}=model.rxns(s2_need);
               k2=k2+1;
            end
        end
    end
end


%% 
isparasit = cellfun( @(X) X ==1 , Species1_parasit(:,2), 'UniformOutput', true )
numParasit = sum(isparasit);
[all,~]=size(Species1_parasit);
rate_of_Parasit_species1=numParasit /all1;

isparasit = cellfun( @(X) X ==1 , Species2_parasit(:,2), 'UniformOutput', true )
numParasit = sum(isparasit);
[all,~]=size(Species2_parasit);
rate_of_Parasit_species2=numParasit /all2;






end