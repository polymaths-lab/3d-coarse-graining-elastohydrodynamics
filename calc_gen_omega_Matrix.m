%Relates rate of change of generator r to angular velocity. (Eqn 2.18)
function D = calc_gen_omega_Matrix(r)
rmag = norm(r);
cross_term = sinc(rmag/pi)^2 * [0 -r(3) r(2); r(3) 0 -r(1); -r(2) r(1) 0];
r1 = r(1); r2 = r(2); r3 = r(3);
term3 = -(rmag-sin(rmag)*cos(rmag))/(rmag^3) * [r2^2+r3^2, -r1*r2, -r1*r3; -r1*r2, r1^2+r3^2, -r2*r3; -r1*r3, -r2*r3, r1^2+r2^2];
term3(isnan(term3)) = 0;
D = 2*eye(3) + 2*cross_term + 2*term3;
end