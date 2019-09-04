function Z_new = change_guess(A,Z,Q)
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
    L = A.dim(1)*unit;
    W = A.dim(2)*unit;
    H = A.dim(3)*unit; 
elseif A.shape == 'spiked_shape'
    L = A.dim(1)*unit;
    W = A.dim(3)*unit;
    H = A.dim(5)*unit;
elseif A.shape == 'curved_shape'
    w = 0;
    L = A.dim(3)*unit;
    H = 0;
elseif A.shape == 'spiked_ended'
    L = A.dim(1)*unit;
    W = A.dim(3)*unit;
    H = A.dim(5)*unit;
elseif A.shape == 'geckod_shape'
    L = A.dim(1)*unit;
    W = A.dim(2)*unit;
    H = A.dim(3)*unit; 
end




% ECP1


R = rand;
a11_x = Z(7);
a11_y = Z(8)+L*R;
a11_z = Z(9);

a12_x = a11_x;
a12_y = a11_y;
a12_z = Z(12);






    % normal force
    if A.shape == 'cuboid_shape'
        p1_n =Z(24);
    elseif  A.shape == 'spiked_shape'
        p1_n =Z(26);
    elseif  A.shape == 'curved_shape'
        p1_n =Z(24);
    elseif  A.shape == 'spiked_ended'
        p1_n =Z(28);
    elseif A.shape == 'geckod_shape'
        p1_n =Z(24);
    end
%     if p1_n == 0
%         a11_z = Z(9)+R*Q(3);
%     else
%         a11_z = 0;
%     end
    
    R = rand;
    p1_n = R*p1_n;

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

    v1_t = v_x - w_z*(a11_x - q_x)+ w_y*(a11_z - q_z);
    v1_o = v_y + w_z*(a11_y - q_y)- w_x*(a11_z - q_z);
    v1_r = w_z;

    sig1 = sqrt(e_t^2*v1_t^2 + e_o^2*v1_o^2 + e_r^2*v1_r^2);
    if sig1 == 0
        p1_t = Z(13);
        p1_o = Z(14);
        p1_r = Z(15);
    else
        p1_t = -e_t^2*mu*p1_n*v1_t/sig1;
        p1_o = -e_o^2*mu*p1_n*v1_o/sig1;
        p1_r = -e_r^2*mu*p1_n*v1_r/sig1;
    end

    

    Z_new = Z;
%     Z_new(1) = v_x;
%     Z_new(5) = 0;
%     Z_new(6) = 0;
    Z_new(7:12) = [a11_x;a11_y;a11_z;a12_x;a12_y;a12_z];
    Z_new(13) = p1_t;
    Z_new(14) = p1_o;
    Z_new(15) = p1_r;
    Z_new(16) = sig1;
    if A.shape == 'cuboid_shape'
        Z_new(24) =0;
    elseif  A.shape == 'spiked_shape'
        Z_new(26) = p1_n+m*g*h*(R);
    elseif  A.shape == 'spiked_ended'
        Z_new(28) = p1_n;
    elseif A.shape == 'geckod_shape'
        Z_new(24) = p1_n;
    elseif A.shape == 'curved_shape'
        Z_new(24) = p1_n;
    end





end