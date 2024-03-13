function [X,params] = get_params_chlamy()
Nstruct = 1; %one structure
Nfil = 2*ones(Nstruct,1); %with two filaments per structure

N = 13; %with N segments per filament
n = 1; %with n spheres per segment

a = 1/(2*N*n);

R = 0.2; %body of radius 0.3 L
ra = 4*a; %body spheres have radius 4x filament spheres here, just to help code run fast for quick example.
%Actual use case set ra to a smaller value to get a nicer looking body. Use
%MEX function for calc_Mh to speed up mobility matrix calculations.

% -------
%This code uses a very simple expression for placing spheres on a larger
%sphere. You may want to let spheres evolve under a repulsive potential for
%long time.
number_of_layers = R*pi/(2*ra);
number_of_layers = ceil(number_of_layers)-1;
theta = 0;
c=0;
for layer = 1:number_of_layers
    Nlayer = pi*R*sin(theta)/(ra);
    Nlayer = floor(Nlayer);
    
    phi = 0;
    for i = 1:Nlayer
        c=c+1;
        b(c,:) = R*[sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
        phi = phi+2*pi/Nlayer;
        
    end
    theta = theta + pi/number_of_layers;
end

b(c+1,:) = [0 0 R]; %add sphere perfectly on top
b(c+2,:) = [0 0 -R]; %and on bottom
%---------

Nbody = size(b,1);
params.Nbody = Nbody;

%Set size of spheres in array 'a'
sn = 0;
for i = 1:Nstruct
    for j = 1:Nbody(i)
        sn = sn + 1;
        a(sn) =  ra;
    end
    
    for j = 1:Nfil(i)
        for k = 1:N*n
            sn = sn + 1;
            a(sn) = 1/(2*N*n);
        end
    end
end


%Nbody and Nfil are arrays that are Nstruct long. Nbody(i) and Nfil(i)
%contain the number of body spheres and number of filaments in struct i.
%This checks that they are the correct length.
if length(Nfil) > Nstruct
    Nfil = Nfil(1);
end
if length(Nbody) > Nstruct
    Nbody = Nbody(1);
end


%Set b for body and filaments of all structures
%Nstruct is 1 here, so outer loop not needed. This tells the code where to
%put the base of filament j relative to x0.
an = 0.2*[-1 1];
for i = 1:Nstruct
    for j = 1:Nfil(i)
        b(Nbody(i)+j,:,i) = [(R+a(1))*sin(an(j)) 0 (R+a(1))*cos(an(j))];
    end
end

params.Nstruct = Nstruct;
params.N = N;
params.n = n;
params.Nbody = Nbody;
params.Nfil = Nfil;
params.b = b;
params.a = a';

X(1:3) = [0 0 0]';
X(4:6) = [0 0 0]';
X = X(:);

%Set how the flagella exit the chlamydomonas body
for fil = 1:Nfil
    u = [0 1 0];
    for i = 1:N
        theta = an(fil);
        u = [0 1 0];
        u = u / norm(u);
        r = (theta/2) * u;
        X(3+3+3*N*(fil-1)+3*(i-1)+1:3+3+3*N*(fil-1)+3*i) = r';
    end
end


end