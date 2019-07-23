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

    
Z_new = Z;

R = rand;    
Z_new(44) = R*Z(44);

Z_new(52) = R*Z(52);

Z_new(60) = R*Z(60);
end