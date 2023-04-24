%Calculates matrix that encodes the force and torque balance of entire
%structure and sub-filaments. Eqn 9, 10
function [A,B,C] = calc_Mf(X, Nbody, Nfil, N,b,n)
%Find the matrix that relates forces and torques on every sphere within a
%particular filament to the elastic forces and torques within the filament
%e.g. Mf * [F; T] = [0, 0, 0, \kappa11, \kappa12, \kappa13, ...]

%Get quaternion elements
q0 = X(3+1:4:end-3);
q1 = X(3+2:4:end-2);
q2 = X(3+3:4:end-1);
q3 = X(3+4:4:end-0);

%Get director basis
d1 = [q0.*q0 + q1.*q1 - q2.*q2 - q3.*q3, 2*(q2.*q1 + q0.*q3), 2*(q3.*q1-q0.*q2)];
d2 = [2*(q1.*q2-q0.*q3), q0.*q0-q1.*q1+q2.*q2-q3.*q3, 2*(q3.*q2+q0.*q1)];
d3 = [2*(q1.*q3+q0.*q2), 2*(q2.*q3-q0.*q1), q0.*q0-q1.*q1-q2.*q2+q3.*q3];

%Total force = 0
A = zeros(3, 3*(Nbody+Nfil*N*n));
% A = repmat(eye(3),Nbody+Nfil*N*n);
for i = 1:Nbody+(Nfil*N*n)  
    A(1:3, 3*(i-1)+1:3*i) = eye(3);
end

x0 = X(1:3); %body of struct
X3 = calc_sphere_centres_full(X, Nbody, Nfil, N,b,n);

B = zeros(3*N*Nfil+3, 3*(Nbody+(Nfil*N*n)));
C = zeros(3*N*Nfil+3, 3*(Nbody+(Nfil*N*n)));

 %total torque on entire structure about x0 = 0
for i = 1:Nbody+(Nfil*N*n)
    xi = X3(3*(i-1)+1:3*i);
    dx = xi - x0;
    dxcross = [0 -dx(3) dx(2);dx(3) 0 -dx(1); -dx(2) dx(1) 0];
    Bji = dxcross;
    
    B(3*(1-1)+1:3*1, 3*(i-1)+1:3*i) = Bji;
    C(3*(1-1)+1:3*1, 3*(i-1)+1:3*i) = eye(3);
end

%subfilament torques
for fil = 1:Nfil
    for j = 1:N
        %sphere number
        snj = Nbody + (fil-1)*N*n+(j-1)*n+1;

        %eqn number
        eqn = 1+(fil-1)*N + j;
        
        d1j = d1(1+(fil-1)*N+j,:)';
        d2j = d2(1+(fil-1)*N+j,:)';
        d3j = d3(1+(fil-1)*N+j,:)';        

        xj = X3(3*(snj-1)+1:3*snj) - d3j*1/(2*N*n);

        for i = (j-1)*n+1:N*n
            sni = Nbody + (fil-1)*N*n+i;
            xi = X3(3*(sni-1)+1:3*sni);
            dx = (xi-xj);
            dxcross = [0 -dx(3) dx(2);dx(3) 0 -dx(1); -dx(2) dx(1) 0];
            Bji = dxcross;
            B(3*(eqn-1)+1:3*eqn, 3*(sni-1)+1:3*sni) = Bji;
            C(3*(eqn-1)+1:3*eqn, 3*(sni-1)+1:3*sni) = eye(3);
        end              
    end
end
end