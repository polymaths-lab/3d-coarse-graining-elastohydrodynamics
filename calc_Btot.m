function Btot = calc_Btot(t,X, params)
Xq = calc_Xq(X,params);
N = params.N;
n = params.n;
Nstruct = params.Nstruct;
Nfil = params.Nfil;
a = params.a;
b = params.b;
Nbody = params.Nbody;
Btot = [];
s=0;
for i = 1:Nstruct
    M(i) = 3+4+4*N*Nfil(i); %size of X for this struct    
    Xstruct = Xq(s + 1:s + M(i));

    B = calc_RHS(t,Xq,params);

    Btot = cat(1,Btot, B);
    s = s + M(i);
end

%
s=0;
for i = 1:Nstruct
    Neqns = 3+3+3*N*Nfil(i);
    
    for fil = 1:Nfil(i)
        j = 1;
        snj = Nbody(i) + (fil-1)*N + j;
        eqns = s+6+(snj-Nbody(i)-1)*3+1:s+6+(snj-Nbody(i))*3; %out of struct i eqns

        Btot(eqns) = 0;
    end
    s = s + 3+3+3*N*Nfil(i);
end

end
    
    
    