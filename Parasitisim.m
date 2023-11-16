function [Isinhibite] = Parasitisim(v_demand,biomass,v_rate,bio_rate,model)

[m,n]=size(model.S);

Aeq=zeros(m+1,n);
beq=zeros(1,m+1);

Aeq(1:m,1:n)=model.S(1:m,1:n);

Aeq(m+1,v_demand)=1;
beq(1,m+1)=v_rate;


A=zeros(1,n);
A(1,biomass)=-1;
b=-1*(1-0.001)*bio_rate;

f=zeros(1,n);

options = optimoptions('linprog','Algorithm','dual-simplex');

[x,fval,exitflag,output] = linprog(f,A,b,Aeq,beq,model.lb,model.ub,options);

if(exitflag == 1)
    Isinhibite=0;
else
    Isinhibite=1;
end

end