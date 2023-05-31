%Calculates matrix relating rate of change of generator state vector for
%structure to the angular velocity of segments. (Eqn 2.21).
function D = calc_D(Xstruct, params, struct)
    N = params.N;
    Nfil = params.Nfil;
    
    D = zeros(3+3+3*N*Nfil(struct),3+3+3*N*Nfil(struct));
    
    r0 = Xstruct(4:6); 
    D(1:3,1:3) = eye(3);
    D(4:6,4:6) =  calc_gen_omega_Matrix(r0);
    s = 0;
    for i = 1:Nfil(struct)
        M(i) = 3*N;
        for j = 1:N
            r = Xstruct(s+6+3*(j-1)+1:s+6+3*j);            
            D(s+6+3*(j-1)+1:s+6+3*j, s+6+3*(j-1)+1:s+6+3*j) = calc_gen_omega_Matrix(r);
        end
        s = s + M(i);
    end
end
