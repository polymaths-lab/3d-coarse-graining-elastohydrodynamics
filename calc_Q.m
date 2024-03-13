%Calculate Q matrix that converts vector of angular velocities into vector
%of velocities and angular velocities, Eqn 2.13.
function Q = calc_Q(X, N,Nfil,Nbody, b,n)
q0 = X(3+1:4:end-3);
q1 = X(3+2:4:end-2);
q2 = X(3+3:4:end-1);
q3 = X(3+4:4:end-0);

d1 = [q0.*q0 + q1.*q1 - q2.*q2 - q3.*q3, 2*(q2.*q1 + q0.*q3), 2*(q3.*q1-q0.*q2)];
d2 = [2*(q1.*q2-q0.*q3), q0.*q0-q1.*q1+q2.*q2-q3.*q3, 2*(q3.*q2+q0.*q1)];
d3 = [2*(q1.*q3+q0.*q2), 2*(q2.*q3-q0.*q1), q0.*q0-q1.*q1-q2.*q2+q3.*q3];
d1=load("d1.txt");
d2=load("d2.txt");
d3=load("d3.txt");
% d3(3,1) = 400
a = 1/(2*N*n);

Q = zeros(6*(Nbody+(Nfil*N*n)), 6+3*N*Nfil);
Q(1:3, 1:3) = eye(3);

for i = 1:Nbody
    r=b(i,:)*[d1(1,:); d2(1,:); d3(1,:)];
    rcross = [0 -r(3) r(2); r(3) 0 -r(1); -r(2) r(1) 0];
    Q(3*(i-1)+1:3*i, 1:3) = eye(3);
    Q(3*(i-1)+1:3*i, 4:6) =  -rcross;
    
    Q(3*(Nbody+(Nfil*N*n))+3*(i-1)+1:3*(Nbody+(Nfil*N*n))+3*i, 4:6) = eye(3);
end
%
for fil = 1:Nfil
    seg_n = Nbody + (fil-1)*N+1;
    %     sph_n = Nbody + (fil-1)*N*n + 1;
            
    r=b(Nbody+fil,:)*[d1(1,:); d2(1,:); d3(1,:)];
    rcross = [0 -r(3) r(2); r(3) 0 -r(1); -r(2) r(1) 0];
    d3i = d3(1+seg_n-Nbody, :);
    
    d3cross = [0 -d3i(3) d3i(2); d3i(3) 0 -d3i(1); -d3i(2) d3i(1) 0];
    
    for k = 1:n
        sn = Nbody+(fil-1)*N*n+k;
        
        Q(3*(sn-1)+1:3*sn, 1:3) =  eye(3);
        Q(3*(sn-1)+1:3*sn, 4:6) = -rcross; %gets to start of filament
        
        Q(3*(sn-1)+1:3*sn, 3+3+(seg_n-1-Nbody)*3+1:3+3+(seg_n-Nbody)*3) = -(2*k*a-a)*d3cross; %v
        
        Q(3*(Nbody+(Nfil*N*n))+3*(sn-1)+1:3*(Nbody+(Nfil*N*n))+3*sn,3+3+3*(seg_n-Nbody-1)+1:3+3+3*(seg_n-Nbody)) = eye(3); %omega
    end
    for i = 2:N        
        seg_n = Nbody + (fil-1)*N+i;
   
        for k = 1:n
            d3i = d3(1+seg_n-Nbody, :);
            d3cross = [0 -d3i(3) d3i(2); d3i(3) 0 -d3i(1); -d3i(2) d3i(1) 0];
            sn = Nbody+(fil-1)*N*n+(i-1)*n+k;
            Q(3*(sn-1)+1:3*sn, 1:3) = eye(3);
            Q(3*(sn-1)+1:3*sn, 4:6) = -rcross;
            
            %set rows = to rows from previous sphere
            Q(3*(sn-1)+1:3*sn, :) = Q(3*(sn-1-1)+1:3*(sn-1), :);
            
            Q(3*(sn-1)+1:3*sn, 3+3+(seg_n-Nbody-1)*3+1:3+3+(seg_n-Nbody)*3) = -(2*k*a-a)*d3cross;
            
            d3i = d3(1+seg_n-Nbody-1, :);
            d3cross = [0 -d3i(3) d3i(2); d3i(3) 0 -d3i(1); -d3i(2) d3i(1) 0];
            Q(3*(sn-1)+1:3*sn, 3+3+(seg_n-Nbody-1-1)*3+1:3+3+(seg_n-Nbody-1)*3) = -1/N*d3cross;
            Q(3*(Nbody+(Nfil*N*n))+3*(sn-1)+1:3*(Nbody+(Nfil*N*n))+3*sn, 3+3+3*(seg_n-Nbody-1)+1:3+3+3*(seg_n-Nbody)) = eye(3);
        end
    end
end
end