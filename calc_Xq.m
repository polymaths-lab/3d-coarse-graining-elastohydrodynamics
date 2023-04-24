%Input generator state vector X and return quaternion state vector Xq.
function Xq = calc_Xq(X,params)
N = params.N;
Nfil = params.Nfil;
Nstruct = params.Nstruct;

s=0;
sq = 0;
Xq = zeros(7*Nstruct+4*N*sum(Nfil),1);
for i = 1:Nstruct
    M(i) = 6+3*N*Nfil(i); %number of unknowns for ith struct
    Mq(i) = 7+4*N*Nfil(i); %number of unknowns for ith struct
    
    Xq(sq+1:sq+3) = X(s+1:s+3);
    r = X(s+4:s+6);
    rmag = norm(r);
    q = [cos(rmag); sinc(rmag/pi)*r];
    Xq(sq+4:sq+7) = q;
    kk = 0;
    for j = 1:Nfil(i)
        for k = 1:N
            kk = kk + 1;
            r = X(s+6+3*(kk-1)+1:s+6+3*kk);
            rmag = norm(r);
            q = [cos(rmag); sinc(rmag/pi)*r];
            Xq(sq+7+4*(kk-1)+1:sq+7+4*kk) = q;
        end
    end
    s = s + M(i);
    sq = sq + Mq(i);
end