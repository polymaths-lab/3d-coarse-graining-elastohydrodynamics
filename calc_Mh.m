%Calculate RPY matrix, relating forces and torques on all spheres to
%linear and angluar velocities of all spheres. Eqn 2.11, and Appendix.
function Mh = calc_Mh(Ntot,X3f,a)

Mh = zeros(6*Ntot, 6*Ntot);

for i = 1:Ntot
    xi = X3f(3*(i-1)+1:3*i);
    for j = i+1:Ntot
        xj = X3f(3*(j-1)+1:3*j);
        rij = (xi - xj);
        r = norm(rij);
        rhat = rij/r;
        rhat_skew = [0 -rhat(3) rhat(2); ...
            rhat(3) 0 -rhat(1); ...
            -rhat(2) rhat(1) 0]';
        
        if r > a(i) + a(j)
            Mtt = 1/(2*r) * ((1 + (a(i)^2 + a(j)^2)/(3*r^2))*eye(3) + (1-(a(i)^2+a(j)^2)/(r^2))*rhat*rhat');
            Mtr = 1/(2*r^2)*rhat_skew;
            Mrt = Mtr;
            Mrr = 1/(4*r^3)*(3*rhat*rhat' - eye(3));
        else
            Mtt=2/(3*a(i)*a(j))*(((16*r^3*(a(i)+a(j))-((a(i)-a(j))^2+3*r^2)^2)/(32*r^3))*eye(3)+ ...
                3*((a(i)-a(j))^2-r^2)^2/(32*r^3)*rhat*rhat');

            temp=(5*r^6-27*r^4*(a(i)^2+a(j)^2)+32*r^3*(a(i)^3+a(j)^3)-9*r^2*(a(i)^2-a(j)^2)^2-(a(i)-a(j))^4*(a(i)^2+4*a(i)*a(j)^2))/(64*r^3);
            temp2=(3*((a(i)-a(j))^2-r^2)^2*(a(i)^2+4*a(i)*a(j)+a(j)^2-r^2))/(64*r^3);
            Mrr=1/(2*a(i)^3*a(j)^3)*(temp*eye(3)+temp2*rhat*rhat');
            Mtr = 1/(4*a(i)^3*a(j))*(a(i)-a(j)+r)^2*(a(j)^2+2*a(j)*(a(i)+r)-3*(a(i)-r)^2)/(8*r^2)*rhat_skew;
            Mrt = Mtr;
        end

        Mh(3*(i-1)+1:3*i, 3*(j-1)+1:3*j) = Mtt;
        Mh(3*(i-1)+1:3*i, 3*(Ntot)+3*(j-1)+1:3*(Ntot)+3*j) = Mtr;
        Mh(3*(Ntot)+3*(i-1)+1:3*(Ntot)+3*i, 3*(j-1)+1:3*j) = Mrt;
        Mh(3*(Ntot)+3*(i-1)+1:3*(Ntot)+3*i, 3*(Ntot)+3*(j-1)+1:3*(Ntot)+3*j) = Mrr;

        Mh(3*(j-1)+1:3*j, 3*(i-1)+1:3*i) = Mtt';
        Mh(3*(j-1)+1:3*j, 3*(Ntot)+3*(i-1)+1:3*(Ntot)+3*i) = Mtr';
        Mh(3*(Ntot)+3*(j-1)+1:3*(Ntot)+3*j, 3*(i-1)+1:3*i) = Mrt';
        Mh(3*(Ntot)+3*(j-1)+1:3*(Ntot)+3*j, 3*(Ntot)+3*(i-1)+1:3*(Ntot)+3*i) = Mrr';
    end

    Mtt = 2/(3*a(i))*eye(3);
    Mrr = 1/(2*a(i)^3)*eye(3);
    Mtr = zeros(3);
    Mrt = Mtr;
    Mh(3*(i-1)+1:3*i, 3*(i-1)+1:3*i) = Mtt;
    Mh(3*(i-1)+1:3*i, 3*(Ntot)+3*(i-1)+1:3*(Ntot)+3*i) = Mtr;
    Mh(3*(Ntot)+3*(i-1)+1:3*(Ntot)+3*i, 3*(i-1)+1:3*i) = Mrt;
    Mh(3*(Ntot)+3*(i-1)+1:3*(Ntot)+3*i, 3*(Ntot)+3*(i-1)+1:3*(Ntot)+3*i) = Mrr;
end
Mh = Mh * 1/(4*pi);
end