function [X,params] = get_params_relaxing_filament();
    params.Nstruct = 1; %One structure
    params.Nfil = 1;  %With one filament
    params.Nbody = 0; %No body

    params.N = 20; %20 segments
    params.n = 2;  %2 spheres per segment
    
    params.a = 1/(2*params.N*params.n)*ones(params.N*params.n,1); %sphere radii
    params.b = [0 0 0]; %xb = x1, separation between xb and x1 is zero
    
    X = zeros(6+3*params.N,1);

    X(1:3) = [0 0 0]; %Position of filament base
    X(4:6) = [0 0 0]; %Body frame, irrelevant here because there is no body
    X(8:3:end) = linspace(0, pi/2, params.N); %bend filament in plane
end