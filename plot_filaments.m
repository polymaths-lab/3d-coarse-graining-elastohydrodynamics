%Plots all filaments in system using streamtube function.
function h=plot_filaments(X,params,flag)
N = params.N;
n = params.n;
Nstruct = params.Nstruct;
Nfil = params.Nfil;
a = params.a;
b = params.b;
Nbody = params.Nbody;
s=0;

if flag == 1
    X = calc_Xq(X,params);
end
    
for i = 1:Nstruct
    M(i) = 7+4*N*Nfil(i);
    Xi=X(s+1:s+M(i));
    
    q0 = Xi(3+1:4:end-3);
    q1 = Xi(3+2:4:end-2);
    q2 = Xi(3+3:4:end-1);
    q3 = Xi(3+4:4:end-0);
    
    d1t = [q0.*q0 + q1.*q1 - q2.*q2 - q3.*q3, 2*(q2.*q1 + q0.*q3), 2*(q3.*q1-q0.*q2)];
    d2t = [2*(q1.*q2-q0.*q3), q0.*q0-q1.*q1+q2.*q2-q3.*q3, 2*(q3.*q2+q0.*q1)];
    d3t = [2*(q1.*q3+q0.*q2), 2*(q2.*q3-q0.*q1), q0.*q0-q1.*q1-q2.*q2+q3.*q3];
    
    X3 = calc_sphere_centres_full(Xi, Nbody(i), Nfil(i), N, b(:,:,i), n, params);
    x = X3(1:3:length(X3)-2);
    y = X3(2:3:length(X3)-1);
    z = X3(3:3:length(X3)-0);
    for j = 1:Nfil(i)
        X3j = X3(Nbody(i)+1+3*N*n*(j-1):Nbody(i)+3*N*n*j);
        d3 = 1/(2*N*n)*d3t(2+1*N*(j-1):1+1*N*j,:);
        d2 = 1/(2*N*n)*d2t(2+1*N*(j-1):1+1*N*j,:);
        d1 = 1/(2*N*n)*d1t(2+1*N*(j-1):1+1*N*j,:);
        verts = {[[x(1+Nbody(i)+1*N*n*(j-1)) - d3(1,1);x(1+Nbody(i)+1*N*n*(j-1):Nbody(i)+1*N*n*j); x(Nbody(i)+1*N*n*j) + d3(end,1)] ...
            [y(1+Nbody(i)+1*N*n*(j-1))-d3(1,2);y(1+Nbody(i)+1*N*n*(j-1):Nbody(i)+1*N*n*j);y(Nbody(i)+1*N*n*j)+d3(end,2)] ...
            [z(1+Nbody(i)+1*N*n*(j-1))-d3(1,3);z(1+Nbody(i)+1*N*n*(j-1):Nbody(i)+1*N*n*j);z(Nbody(i)+1*N*n*j)+d3(end,3)]]};

        %Plot tube representing filament of diameter 1/(Nn)
        h(i)=streamtube(verts,1/(N*n));
        hold on
        
        h(i).EdgeColor = 'none';
        h(i).FaceColor = [0.9290 0.6940 0.1250];
        h(i).AmbientStrength = 0.5;
        
        as = 1*[0 1 0 -1 0.707 0.707 -0.707 -0.707, 0];
        bs = 1*[1 0 -1 0 0.707 -0.707 0.707 -0.707,0];
        
        xx = [x(1+Nbody(i)+1*N*n*(j-1))-d3(1,1);x(1+1+Nbody(i)+1*N*n*(j-1):Nbody(i)+1*N*n*j-1); x(Nbody(i)+1*N*n*j) + d3(end,1)];
        yy = [y(1+Nbody(i)+1*N*n*(j-1))-d3(1,2);y(1+1+Nbody(i)+1*N*n*(j-1):Nbody(i)+1*N*n*j-1); y(Nbody(i)+1*N*n*j) + d3(end,2)];
        zz = [z(1+Nbody(i)+1*N*n*(j-1))-d3(1,3);z(1+1+Nbody(i)+1*N*n*(j-1):Nbody(i)+1*N*n*j-1); z(Nbody(i)+1*N*n*j) + d3(end,3)];

        %Plot lines on surface of tube to show twist
        for iii = 1:8
            a1 = as(iii);
            b1 = bs(iii);
            li = [xx + a1 * d1(:,1) + b1*d2(:,1),yy + a1 * d1(:,2) + b1*d2(:,2), zz + a1 * d1(:,3) + b1*d2(:,3)];
            plot3(li(:,1),li(:,2), li(:,3),'-',  'LineWidth', 0.2,'Color', 'k');
        end
    end
    
    s = s + M(i);
end
end