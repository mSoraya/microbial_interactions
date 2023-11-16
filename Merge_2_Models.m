function [model] = Merge_2_Models(model1,model2)

[m1,n1]=size(model1.S);
[m2,n2]=size(model2.S);

EX_1=find(contains(model1.rxns,'EX_'));
EX_2=find(contains(model2.rxns,'EX_'));

[i1,~]=find(model1.S(:,EX_1)~=0);
[i2,~]=find(model2.S(:,EX_2)~=0);

mets_1=model1.mets(i1);
mets_2=model2.mets(i2);

shared_mets=intersect(mets_1,mets_2);

clear mets_1 mets_2;

model.S(1:m1,1:n1)=model1.S(1:m1,1:n1);
model.rxns(1:n1,1)=strcat(model1.rxns(1:n1),'_species1');
model.mets(1:m1,1)=strcat(model1.mets(1:m1),'_species1');
model.lb(1:n1,1)=model1.lb(1:n1);
model.ub(1:n1,1)=model1.ub(1:n1);
model.minFlux(1:n1,1)=model1.minFlux(1:n1);
model.maxFlux(1:n1,1)=model1.maxFlux(1:n1);

model.S(m1+1:m1+m2,n1+1:n2+n1)=model2.S(1:m2,1:n2);
model.rxns(n1+1:n1+n2,1)=strcat(model2.rxns(1:n2),'_species2');
model.mets(m1+1:m1+m2,1)=strcat(model2.mets(1:m2),'_species2');
model.lb(n1+1:n1+n2)=model2.lb(1:n2);
model.ub(n1+1:n1+n2)=model2.ub(1:n2);
model.minFlux(n1+1:n1+n2)=model2.minFlux(1:n2);
model.maxFlux(n1+1:n1+n2)=model2.maxFlux(1:n2);


k=n1+n2+1;
for i=1:length(shared_mets)
    clear idx1 idx2 rxn1 rxn2;
    model.mets(m1+m2+i,1)=strcat(extractBefore(shared_mets(i),'['),'_[lu]');

    idx1=find(strcmp(model1.mets,shared_mets(i)));
    idx2=find(strcmp(model2.mets,shared_mets(i)));

    temp=find(model1.S(idx1,:)~=0);
    rxn1=intersect(EX_1,temp);
    clear temp;
    temp=find(model2.S(idx2,:)~=0);
    rxn2=intersect(EX_2,temp);

    model.rxns(rxn1,1)=strcat( 'Re_',model.rxns(rxn1,1));
    model.rxns(rxn2+n1)=strcat( 'Re_',model.rxns(rxn2+n1));


    model.shared_rxn(i,1:length(rxn1)+length(rxn2))=union(rxn1,rxn2+n1);
    model.shared_met(i,1)=m1+m2+i;

    for j=1:length(rxn1)
        model.S(m1+m2+i,rxn1(j))=model1.S(idx1,rxn1(j))*-1;
    end

    for j=1:length(rxn2)
        model.S(m1+m2+i,rxn2(j)+n1)=model2.S(idx2,rxn2(j))*-1;
    end

    if(length(rxn1)+length(rxn2)<3)
        if (sign(model.S(idx1,rxn1))==sign(model.S(idx2+m1,rxn2+n1)))
            model.S(m1+m2+i,k)=sign(model.S(idx1,rxn1))*1;

            if(sign(model.S(idx1,rxn1))==1)
                model.rxns(k,1)=strcat(extractBefore(shared_mets(i),'['),'_f_EX_[lu]');
            else
                model.rxns(k,1)=strcat(extractBefore(shared_mets(i),'['),'_b_EX_[lu]');
            end
            model.lb(k,1)=0;
            model.ub(k,1)=1000;
            model.minFlux(k,1)=0;
            model.maxFlux(k,1)=1000;

            k=k+1;
        else
            model.S(m1+m2+i,k)=1;
            model.rxns(k,1)=strcat(extractBefore(shared_mets(i),'['),'_b_EX_[lu]');
            model.lb(k,1)=0;
            model.ub(k,1)=1000;
            model.minFlux(k,1)=0;
            model.maxFlux(k,1)=1000;

            k=k+1;

            model.S(m1+m2+i,k)=-1;
            model.rxns(k,1)=strcat(extractBefore(shared_mets(i),'['),'_f_EX_[lu]');
            model.lb(k,1)=0;
            model.ub(k,1)=1000;
            model.minFlux(k,1)=0;
            model.maxFlux(k,1)=1000;
            k=k+1;
        end
    else
        model.S(m1+m2+i,k)=1;
        model.rxns(k,1)=strcat(extractBefore(shared_mets(i),'['),'_b_EX_[lu]');
        model.lb(k,1)=0;
        model.ub(k,1)=1000;
        model.minFlux(k,1)=0;
        model.maxFlux(k,1)=1000;

        k=k+1;

        model.S(m1+m2+i,k)=-1;
        model.rxns(k,1)=strcat(extractBefore(shared_mets(i),'['),'_f_EX_[lu]');
        model.lb(k,1)=0;
        model.ub(k,1)=1000;
        model.minFlux(k,1)=0;
        model.maxFlux(k,1)=1000;
        k=k+1;
    end
end


end