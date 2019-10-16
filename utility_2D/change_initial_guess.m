function Z_new = change_initial_guess(A,Z)
unit = A.unit;
unit_mass = A.unit_mass;
e_t = A.ellipsoid(1);
e_o = A.ellipsoid(2);
e_r = A.ellipsoid(3)*unit;

m = A.mass*unit_mass;
g = A.gravity*unit;
h = A.h;
mu = A.cof;

R = rand;    
Z_new = R*Z;

R = rand;    
Z_new(48) = R*Z(48);

Z_new(49) = R*Z(49);

Z_new(50) = R*Z(50);

Z_new(51) = R*Z(51);
end