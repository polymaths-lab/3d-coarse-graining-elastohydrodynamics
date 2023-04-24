function h=plot_bodies(X,params,flag)
N = params.N;
n = params.n;
Nstruct = params.Nstruct;
Nfil = params.Nfil;
a = params.a;
b = params.b;
Nbody = params.Nbody;
Xq = X;
h = [];
if flag == 1
    Xq = calc_Xq(X,params);
end

%Get position of all spheres in system
X3 = [];
s=0;
for i = 1:Nstruct
    M(i) = 7+4*N*Nfil(i);
    Xqq=Xq(s+1:s+M(i));
    X3i=calc_sphere_centres_full(Xq(s+1:s+M(i)), Nbody(i), Nfil(i), N, b(:,:,i), n);
    X3 = cat(1,X3,X3i(1:3*Nbody(i)));
    s = s + M(i);
end
[A,B,C] = sphere(10);

%Plot just spheres that are part of bodies
c = 1;
j = 0;
k = 0;
for i = 1:length(X3)/3
    k = k + 1;
    if k > Nbody(c)
        j = j + N*n*Nfil(c)+Nbody(c);
        c = c + 1;
        k=1;
    end
    h(i) = surf(X3((i-1)*3+1)+a(j+k)*A,X3((i-1)*3+2)+a(j+k)*B,X3((i-1)*3+3)+a(j+k)*C);
    hold on
    h(i).FaceColor = [0.8500 0.3250 0.0980];
 
end
end