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
    params.b(:,:,2) = [0 0 0];
    
    X = zeros(6+3*params.N,1);

    X(1:3) = [0 0 0]; %Position of filament base in lab frame
    X(4:6) = [0 0 0]; %Body frame, irrelevant here because there is no body
    X(8:3:end) = linspace(0, pi/2, params.N); %bend filament in plane

%     X = [X;X];
%     X(1:3) = [-1 -1 -1];
end