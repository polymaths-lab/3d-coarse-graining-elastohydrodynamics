function [X,params] = get_params_relaxing_filament();
    params.Nstruct = 1; %One structure
    params.Nfil = ones(params.Nstruct,1);  %With one filament
    params.Nbody = zeros(params.Nstruct,1); %No body

    params.N = 20; %number of segments constructing filaments
    params.n = 1;  %number of spheres constructing segments
    
    Ntot = sum(params.Nfil*params.N*params.n)+sum(params.Nbody);
    params.a = 1/(2*params.N*params.n)*ones(Ntot,1); %sphere radii
    
    %Position of filament base in body frame. 
    %xb = x1, separation between xb (body frame) and x1 (segment 1 frame) is zero
    params.b = [0 0 0];

    
    X = zeros(6+3*params.N,1);

    %6 + 3N state vector. X(1:3) contains position of base of filament
    %X(4:6) contains body generator. Generally dX(4:6)/dt = d(X(7:9))/dt.
    %X(7:end) contains the generators/exponential map for all segments.
    X(1:3) = [0 0 0]; %Position of filament base in lab frame
    X(4:6) = [0 0 0]; %Body frame, irrelevant here because there is no body

    X(7:3:end) = 0.0 %first component of generators
    X(8:3:end) = linspace(0, pi/2, params.N); %second component of generators, here we bend filament in plane
    X(9:3:end) = 0.0 %third component of generators

end
