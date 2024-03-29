%Calculates RHS for each structure, using bending moment relation Eqn 2.8.
%Add intrinsic curvatures for swimming or non-straight reference states
function B = calc_RHS(t,Xq, Nfil, N, i)
q0 = Xq(3+1:4:end-3);
q1 = Xq(3+2:4:end-2);
q2 = Xq(3+3:4:end-1);
q3 = Xq(3+4:4:end-0);

d1 = [q0.*q0 + q1.*q1 - q2.*q2 - q3.*q3, 2*(q2.*q1 + q0.*q3), 2*(q3.*q1-q0.*q2)];
d2 = [2*(q1.*q2-q0.*q3), q0.*q0-q1.*q1+q2.*q2-q3.*q3, 2*(q3.*q2+q0.*q1)];
d3 = [2*(q1.*q3+q0.*q2), 2*(q2.*q3-q0.*q1), q0.*q0-q1.*q1-q2.*q2+q3.*q3];

s = 1/N:1/N:1-1/N;

B = zeros(3+3+3*N*Nfil,1);


for fil = 1:Nfil
    d1f = d1(1+(fil-1)*(N)+1:1+fil*(N),:);
    d2f = d2(1+(fil-1)*(N)+1:1+fil*(N),:);
    d3f = d3(1+(fil-1)*(N)+1:1+fil*(N),:);

    d1s=d1f(2:end,:)-d1f(1:end-1,:);
    d2s=d2f(2:end,:)-d2f(1:end-1,:);
    d3s=d3f(2:end,:)-d3f(1:end-1,:);
    s = (1/N:1/N:1-1/N)';

    %Add travelling wave of curvature in d2 direction
    prefk1 = 0;
    prefk2 = 3*sin(2*pi*s-2*pi*t);
    prefk3 = 0;

    %Calculate curvatures in director basis
    k1 =   (N*dot(d2s,d3f(2:end,:),2)-prefk1);
    k2 =   (N*dot(d3s,d1f(2:end,:),2)-prefk2);
    k3 =   (N*dot(d1s,d2f(2:end,:),2)-prefk3);

    %Convert to lab frame
    kx = k1.*d1f(2:end,1)+k2.*d2f(2:end,1)+k3.*d3f(2:end,1);
    ky = k1.*d1f(2:end,2)+k2.*d2f(2:end,2)+k3.*d3f(2:end,2);
    kz = k1.*d1f(2:end,3)+k2.*d2f(2:end,3)+k3.*d3f(2:end,3);

    eqns_s = 6+((fil-1)*N+2-1)*3+1:6+((fil-1)*N+2)*3;
    eqns_f = 6+((fil-1)*N+N-1)*3+1:6+((fil-1)*N+N)*3;

    %Place into RHS
    B(eqns_s(1):3:eqns_f(1)) = - kx;
    B(eqns_s(2):3:eqns_f(2)) = - ky;
    B(eqns_s(3):3:eqns_f(3)) = - kz;
end
end
