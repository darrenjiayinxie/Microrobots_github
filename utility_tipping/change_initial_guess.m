function Z_new = change_initial_guess(A,Z,Q)
unit = A.unit;
unit_mass = A.unit_mass;
e_t = A.ellipsoid(1);
e_o = A.ellipsoid(2);
e_r = A.ellipsoid(3)*unit;

m = A.mass*unit_mass;
g = A.gravity*unit;
h = A.h;
mu = A.cof;

if A.shape == 'cuboid_shape'
    % ECPs
    a1_x = Z(7);
    a1_y = Z(8);
    a1_z = Z(9);



    % normal force
    p_n =Z(24);
elseif  A.shape == 'spiked_shape'
        % ECPs
    a1_x = Z(7);
    a1_y = Z(8);
    a1_z = Z(9);

    % normal force
    p_n =Z(26);
elseif  A.shape == 'spiked_ended'
        % ECPs
    a1_x = Z(7);
    a1_y = Z(8);
    a1_z = Z(9);

    % normal force
    p_n =Z(28);
end


R = rand;
p_n = R*p_n;

nu = Z(1:6); 
v_x = nu(1);
v_y = nu(2);
v_z = nu(3);
w_x = nu(4);
w_y = nu(5);
w_z = nu(6);

q_x = Q(1);
q_y = Q(2);
q_z = Q(3);

v_t = v_x - w_z*(a1_x - q_x)+ w_y*(a1_z - q_z);
v_o = v_y + w_z*(a1_y - q_y)- w_x*(a1_z - q_z);
v_r = w_z;
    
sig = sqrt(e_t^2*v_t^2 + e_o^2*v_o^2 + e_r^2*v_r^2);
if sig == 0
    p_t = Z(13);
    p_o = Z(14);
    p_r = Z(15);
else
    p_t = -e_t^2*mu*p_n*v_t/sig;
    p_o = -e_o^2*mu*p_n*v_o/sig;
    p_r = -e_r^2*mu*p_n*v_r/sig;
end
    
Z_new = Z;
Z_new(13) = p_t;
Z_new(14) = p_o;
Z_new(15) = p_r;
Z_new(16) = sig;
if A.shape == 'cuboid_shape'
    Z_new(24) = p_n;
elseif  A.shape == 'spiked_shape'
    Z_new(26) = p_n;
elseif  A.shape == 'spiked_ended'
    Z_new(28) = p_n;
end
    



end