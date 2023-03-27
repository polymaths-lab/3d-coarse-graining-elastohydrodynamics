%{
Simulate elastic filaments and solid bodies immersed in a low Reynolds
number fluid. Filaments are comprised of N rigid segments with n spheres
per segment. Bodies are comprised of Nbody spheres.
%}

clear
clf

%Get initial condition and simulation parameters.
[X, params] = get_params_relaxing_filament();

%Plot initial condition
h=plot_filaments(X,params,1,0);
if sum(params.Nbody > 0)
    h2 = plot_bodies(X,params,1);
    for i = 1:length(h2)
        h2(i).EdgeColor = 'none';
    end
end
axis equal
view(3)
camlight

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
    
    sol = ode15s(func, [t t+2], X, ode_ops); %solve for small section of time, after which we will check if any generators are near singularities
    
    t = sol.x(end);
    X = sol.y(:,end);

    %Check if any generators need rescaling
    s=0;
    for i = 1:params.Nstruct
        M(i) = 6+3*params.Nfil(i)*params.N;
        if norm(X(s+4:s+6)) > pi/2
            X(s+4:s+6) = X(s+4:s+6) - pi*X(s+4:s+6)/norm(X(s+4:s+6));
        end
        for j = 1:params.Nfil(i)*params.N
            if norm(X(s+6+(j-1)*3+1:s+6+3*j)) > pi/2
                X(s+6+(j-1)*3+1:s+6+3*j) = X(s+6+(j-1)*3+1:s+6+3*j) - pi*X(s+6+(j-1)*3+1:s+6+3*j)/norm(X(s+6+(j-1)*3+1:s+6+3*j));
            end
        end
        s = s+M(i);
    end

    s = sol.x(1);
    f = sol.x(end);

    ts = find(tps>=s,1);
    tf = find(tps<=f,1,'last');
    if ts > tf
        continue
    end

    %Deval on time points in tps array
    traj(ts:tf, :) = deval(sol, tps(ts:tf))';

    %Display current time
    disp(t)
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

        X3 = calc_sphere_centres_full(Xstruct, params.Nbody,params.Nfil,params.N,params.b,params.n);
        X3f = cat(1,X3f, X3);
    end
end
off = [-0.3 0.3];
x_limits = [min(X3f(1:3:end)) max(X3f(1:3:end))] + off;
y_limits = [min(X3f(2:3:end)) max(X3f(2:3:end))] + off;
z_limits = [min(X3f(3:3:end)) max(X3f(3:3:end))] + off;


%Plot solution over full time range
step = 2;
c = 0;

for ti = 1:length(tps)
    clf
    X = traj(ti,:)';

    h=plot_filaments(X,params,1,t);
    if sum(params.Nbody > 0)
        h2 = plot_bodies(X,params,1);
    end

    title(sprintf('Time = %.2f', tps(ti)));
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