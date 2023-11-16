function [IsCompete] = Competition(v1,v2 ,model)

[m,n]=size(model.S);

Exs=find(contains(model.rxns,'EX_[lu]') | startsWith(model.rxns,'EX_'));

intcon=[n+1:n+length(Exs)];

n_v=length(Exs)+n;
n_eq=m+2;

Aeq=zeros(n_eq,n_v);

Aeq(1:m,1:n)=model.S;
Aeq(m+1,v1)=1;
Aeq(m+2,v2)=1;

beq=zeros(1,n_eq);

lb=zeros(1,n_v);
lb(1,1:n)=model.lb;
ub(1,1:n)=model.ub;
ub(1,n+1:n_v)=1;

ub(1,v1)=0;
ub(1,v2)=0;

A=zeros(length(Exs)+1,n_v);
b=zeros(1,length(Exs)+1);

for i=1: length(Exs)
    A(i,Exs(i,1))=-1;
    A(i,n+i)=0.01;
end

A(end,n+1:n+length(Exs))=-1;
b(1,end)=-1;
f=[];

options = optimoptions('intlinprog','Heuristics','rss');

[x,fval,exitflag,output] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub,options);

if(exitflag==2 | exitflag==1)
    IsCompete=0;
else

    IsCompete=1;

end

end