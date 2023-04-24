%Simulate elastic filaments and solid bodies immersed in a low Reynolds
%number fluid. Filaments are comprised of N rigid segments with n spheres
%per segment. Bodies are comprised of Nbody spheres.


%Get initial condition and simulation parameters.
[X, params] = get_params_relaxing_filament();

%Plot initial condition
h = plot_filaments(X,params,1);
if sum(params.Nbody > 0)
    h2 = plot_bodies(X,params,1);
    for i = 1:length(h2)
        h2(i).EdgeColor = 'none';
    end
end
axis equal
view(3)
axis tight

%%

%Set initial time and time simulation time
t = 0; %initial time
T = 10; %simulation time

%Set dimensionless stiffness parameter S = L^4 \mu / (E_b T)
S = 4;

tic
%Solve from t to T
tps = linspace(0,T,T*10); %timepoints to save solution at
traj = zeros(length(tps),length(X)); %save solution in traj matrix
while t < T
    func = @(t, X) calc_Btot(t, X, params); %set RHS function.
    mass_func = @(t,X) mass_m(t, X, params, S); %set mass matrix function.
    ode_ops = odeset('Mass', mass_func); %tell MATLAB solver to use mass matrix

    sol = ode15s(func, [t t+2], X, ode_ops); %solve for small section of time, after which we will check if any generators are near singularities, could also use event condition to check

    t = sol.x(end);
    X = sol.y(:,end);

    %Check if any generators need rescaling
    s=0;
    for i = 1:params.Nstruct %loop over all structures in system
        M(i) = 6+3*params.Nfil(i)*params.N;
        if norm(X(s+4:s+6)) > pi/2 %check body generator
            X(s+4:s+6) = X(s+4:s+6) - pi*X(s+4:s+6)/norm(X(s+4:s+6));
        end
        for j = 1:params.Nfil(i)*params.N %check all filament generators in structure i
            if norm(X(s+6+(j-1)*3+1:s+6+3*j)) > pi/2
                X(s+6+(j-1)*3+1:s+6+3*j) = X(s+6+(j-1)*3+1:s+6+3*j) - pi*X(s+6+(j-1)*3+1:s+6+3*j)/norm(X(s+6+(j-1)*3+1:s+6+3*j));
            end
        end
        s = s + M(i);
    end

    %Get start and end time that have been solved for
    s = sol.x(1);
    f = sol.x(end);

    %Find indices in tps array that we can deval
    ts = find(tps>=s,1);
    tf = find(tps<=f,1,'last');

    if ts > tf
        continue
    end

    %Deval on time points in tps array
    traj(ts:tf, :) = deval(sol, tps(ts:tf))';

    %Display current time
    fprintf('t = %.2f, %.1f%% completed \n',t,t/T*100)
end
toc
%%
%Find x, y and z limits of solution for plotting
for ti = 1:round(length(tps)/10):length(tps)
    X = traj(ti,:)';
    Xq = calc_Xq(X, params);

    X3f = [];
    s=0;
    for i = 1:params.Nstruct
        M(i) = 3+4+4*params.N*params.Nfil(i);
        Xstruct = Xq(s + 1:s + M(i));
        s = s + M(i);

        X3 = calc_sphere_centres_full(Xstruct, params.Nbody(i),params.Nfil(i),params.N,params.b(:,:,i),params.n);
        X3f = cat(1,X3f, X3);
    end
end

%Add some extra to axis limits
off = [-0.3 0.3];
x_limits = [min(X3f(1:3:end)) max(X3f(1:3:end))] + off;
y_limits = [min(X3f(2:3:end)) max(X3f(2:3:end))] + off;
z_limits = [min(X3f(3:3:end)) max(X3f(3:3:end))] + off;


%Plot solution over full time range
figure; set(gcf,'color','white')
c = 0;
clear F
for ti = 1:length(tps)
    clf
    X = traj(ti,:)';
    
    h = plot_filaments(X,params,1);
    h2 = plot_bodies(X,params,1);

    title(sprintf('t = %.2f', tps(ti)));
    axis equal;
    axis([x_limits y_limits z_limits]);
    view(3)

    set(gcf,'color','white')
    drawnow;

    %Store frames to create video:
    %c = c + 1;
    %set(gcf,'color','white');
    %F(c) = getframe(gcf);
end