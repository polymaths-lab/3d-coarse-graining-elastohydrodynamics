function M = mass_m(t, X, params, Sp)
%{
Calculates entire mass matrix that multiplies the time derivative of the 
state vector in Eqn. 22. 
%}
N = params.N;
Nbody = params.Nbody;
Nfil = params.Nfil;
b = params.b;
n = params.n;
Nstruct = params.Nstruct;
a = params.a;

%Use exponential map to convert generators to quaternions
Xq = calc_Xq(X,params);

Ntot = sum(Nbody) + sum(Nfil)*N*n;

%Series of counting variables for indexing matrices
s=0;
%Indexing Mf matrix
rf = 0;
cf = 0;
%Indexing Q matrix
rq =0;
cq =0;
%Indexing D matrix
sd = 0;

Btot = [];
X3f = [];

%loop over all structures in the system and calculate matrices
for i = 1:Nstruct  
    %Size of Xq for this structure
    M(i) = 3+4+4*N*Nfil(i); 
    %Get Xq of this strucuture
    Xstruct = Xq(s + 1:s + M(i)); 
    %Number of sphere in this structure
    Ntot_struct = Nbody(i)+Nfil(i)*N*n;

    %Get compoenents of Mf matrix for this structure, and place into
    %general Mftot matrix for the entire system
    [A,B,C] = calc_Mf(Xstruct, Nbody(i), Nfil(i), N, b(:,:,i), n);
    
    Mftot(rf+1:rf+3+3+3*N*Nfil(i), cf+1:cf+3*(Nbody(i)+Nfil(i)*N*n)) = [A; B];
    Mftot(rf+1:rf+3+3+3*N*Nfil(i), 3*Ntot+cf+1:3*Ntot+cf+3*(Nbody(i)+Nfil(i)*N*n)) = [zeros(3, 3*(Nbody(i)+(Nfil(i)*N*n))); C];

    rf = rf + (3+3+3*N*Nfil(i)); %number of rows for this struct    
    cf = cf + 3*(Nbody(i)+Nfil(i)*N*n); %number of columns for this struct
    
    %Calculate Q matrix for this structure, and place into general Qtot
    %matrix for the entire system. Q matrix converts angular velocities of
    %segments into velocities and angular velocities of all spheres.
    Q = calc_Q(Xstruct, N,Nfil(i),Nbody(i), b(:,:,i), n);        
    Q1 = Q(1:3*Ntot_struct,:);
    Q2 = Q(3*Ntot_struct+1:end,:);
    Qtot(rq+1:rq+size(Q1,1), cq+1:cq+size(Q1,2)) = Q(1:3*Ntot_struct,:);
    Qtot(3*Ntot+rq+1:3*Ntot+rq+size(Q1,1), cq+1:cq+size(Q1,2)) = Q(3*Ntot_struct+1:end,:);
    
    %Update index counters for Qtot matrix.
    rq = rq + size(Q1,1);
    cq = cq + size(Q1,2);

    %Calculate D matrix for this structure, and then place into Dtot matrix
    %for entire system. D matrix converts time derivatives of segment generators
    %into angulalar velocities of segments.
    MD(i) = 3+3+3*N*Nfil(i);
    D = calc_D(X,params,i);
    Dtot(sd+1:sd+MD(i), sd+1:sd+MD(i)) = D;
    sd = sd + MD(i);
    
    %Calculate lab frame positions of all spheres in structure
    X3 = calc_sphere_centres_full(Xstruct,Nbody(i),Nfil(i),N,b(:,:,i),n);    
    X3f = cat(1,X3f, X3);
    
    %Update counter to find Xq of the next structure.
    s = s + M(i);
end

Mh = calc_Mh(Ntot,X3f,a');
% Mh = calc_Mh_mex(Ntot,X3f,a');

M = Sp^4 * Mftot * inv(Mh) * Qtot * Dtot;

%Rigidly attach the first segment in filaments to the body that the
%filament is attached to. This results in the first segment and the body
%having equal angular velocities, which we enforce by replacing three rows
%of the mass matrix with conditions on the body and segment generators.
%Condition: w_body - w_1 = 2D(r_body)d(r_body)/dt - 2D(r_1)d(r_1)/dt = 0.
%Angular velocities w and generators r, with matrix D converting between
%generator time derivatives and angular velocity.
s=0;
for i = 1:Nstruct
    Neqns = 3+3+3*N*Nfil(i);    
    for fil = 1:Nfil(i)
        j = 1;
        snj = Nbody(i) + (fil-1)*N + j;
        eqns = s+6+(snj-Nbody(i)-1)*3+1:s+6+(snj-Nbody(i))*3; %out of struct i eqns
        h = size(M, 2);
        M(eqns,:) = zeros(3, h);
        
        r0 = X(s+4:s+6);
        r  = X(eqns);
        D0 = calc_gen_omega_Matrix(r0);
        D  = calc_gen_omega_Matrix(r);

        r0mag = norm(r0);
        rmag  = norm(r);
                
        M(eqns, eqns) = D;
        M(eqns, s+4:s+6) = -D0;
    end
    s = s + 3+3+3*N*Nfil(i);
end
end