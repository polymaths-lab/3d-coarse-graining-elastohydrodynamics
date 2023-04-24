%Calculates position of all spheres in structure with quaternion state
%vector X. X contains the body position, body quaternion and quaternions
%for all segments in all filaments.
function X3 = calc_sphere_centres_full(X, Nbody, Nfil, N, b, n,params)
X3 = zeros((Nbody+N(1)*n*Nfil)*3,1);
a = 1/(2*N*n); %sphere separation

q0 = X(3+1:4:end-3);
q1 = X(3+2:4:end-2);
q2 = X(3+3:4:end-1);
q3 = X(3+4:4:end-0);

d1 = [q0.*q0 + q1.*q1 - q2.*q2 - q3.*q3, 2*(q2.*q1 + q0.*q3), 2*(q3.*q1-q0.*q2)];
d2 = [2*(q1.*q2-q0.*q3), q0.*q0-q1.*q1+q2.*q2-q3.*q3, 2*(q3.*q2+q0.*q1)];
d3 = [2*(q1.*q3+q0.*q2), 2*(q2.*q3-q0.*q1), q0.*q0-q1.*q1-q2.*q2+q3.*q3];


x0 = X(1:3); %body of struct

for i = 1:Nbody
    X3(3*(i-1)+1:3*i) = x0 + (b(i,:)*[d1(1,:); d2(1,:); d3(1,:)])'; %position of body sphere i
end

for i = 1:Nfil    
    xfil0 = x0 + (b(Nbody+i,:)*[d1(1,:); d2(1,:); d3(1,:)])';
    xnext = xfil0;

    for j = 1:N
        for k = 1:n
            sn = (j-1)*n+k;

            X3(3*Nbody+(i-1)*3*N*n+(sn-1)*3+1:3*Nbody+(i-1)*3*N*n+3*sn) = xnext + (2*k*a-a)*d3(1+(i-1)*N+j,:)';
        end
        xnext = xnext + d3(1+(i-1)*N(1) + j,:)' * 1/N(1);
    end
end
end