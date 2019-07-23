function Z_new = change_guess(A,Z,Q,r1_y,r2_y,rc_1,rc_2)
unit = A.unit;
unit_mass = A.unit_mass;
e_t = A.ellipsoid(1);
e_o = A.ellipsoid(2);
e_r = A.ellipsoid(3)*unit;

m = A.mass*unit_mass;
g = A.gravity*unit;
h = A.h;
mu = A.cof;

L = A.dim(1)*unit;
W = A.dim(2)*unit;
H = A.dim(3)*unit; 

r1 = rc_1;
r1_z = 0;

r2 = rc_2;
r2_z = 0;


% ECP1


R = rand;
a11_x = Z(7);
a11_y = Z(8)+L*R;
a11_z = Z(9);

a12_x = a11_x;
a12_y = Z(11);
a12_z = Z(12);

% ECP2
a21_x = a11_x;


a21_y = Z(14);
a21_z = Z(15);

a22_x = a21_x;

a22_z = Z(18);

if rand >= 0.5
    a22_y = a21_y;
else
    a22_y = -sqrt(r1^2  - (a22_z - r1_z)^2)+r1_y;
end

% ECP3
a31_x = a11_x;


a31_y = Z(20);
a31_z = Z(21);

a32_x = a31_x;

a32_z = Z(24);

if rand >= 0.5
    a32_y = a31_y;
else
    a32_y = -sqrt(r2^2  - (a32_z - r2_z)^2)+r2_y;
end




% normal force
p1_n =Z(44);
if rand >= 0.5
    R = rand;
    p1_n = R*1e-3;
else
    p1_n = 0;
end


p2_n =Z(52);
if rand >= 0.5
    R = rand;
    p2_n = R*1e-3;
else
    p2_n = 0;
end

p3_n =Z(60);
if rand >= 0.5
    R = rand;
    p3_n = R*1e-3;
else
    p3_n = 0;
end



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

if sig1 <= 1e-8
    p1_t = Z(25);
    p1_o = Z(26);
    p1_r = Z(27);
else
    p1_t = -e_t^2*mu*p1_n*v1_t/sig1;
    p1_o = -e_o^2*mu*p1_n*v1_o/sig1;
    p1_r = -e_r^2*mu*p1_n*v1_r/sig1;
end

v2_t = v_x - w_z*(a21_x - q_x)+ w_y*(a21_z - q_z);
v2_o = v_y + w_z*(a21_y - q_y)- w_x*(a21_z - q_z);
v2_r = w_z;

sig2 = sqrt(e_t^2*v2_t^2 + e_o^2*v2_o^2 + e_r^2*v2_r^2);

if sig2 <= 1e-8
    p2_t = Z(28);
    p2_o = Z(29);
    p2_r = Z(30);
else
    p2_t = -e_t^2*mu*p2_n*v2_t/sig2;
    p2_o = -e_o^2*mu*p2_n*v2_o/sig2;
    p2_r = -e_r^2*mu*p2_n*v2_r/sig2;
end

v3_t = v_x - w_z*(a31_x - q_x)+ w_y*(a31_z - q_z);
v3_o = v_y + w_z*(a31_y - q_y)- w_x*(a31_z - q_z);
v3_r = w_z;

sig3 = sqrt(e_t^2*v3_t^2 + e_o^2*v3_o^2 + e_r^2*v3_r^2);

if sig3 <= 1e-8
    p3_t = Z(31);
    p3_o = Z(32);
    p3_r = Z(33);
else
    p3_t = -e_t^2*mu*p3_n*v3_t/sig3;
    p3_o = -e_o^2*mu*p3_n*v3_o/sig3;
    p3_r = -e_r^2*mu*p3_n*v3_r/sig3;
end


Z_new = Z;
Z_new(7:12) = [a11_x;a11_y;a11_z;a12_x;a12_y;a12_z];
Z_new(13:18) = [a21_x;a21_y;a21_z;a22_x;a22_y;a22_z];
Z_new(19:24) = [a31_x;a31_y;a31_z;a32_x;a32_y;a32_z];


Z_new(25:27) = [p1_t;p1_o;p1_r];
Z_new(28:30) = [p2_t;p2_o;p2_r];
Z_new(31:33) = [p3_t;p3_o;p3_r];

Z_new(34:36) = [sig1;sig2;sig3];
Z_new(37:43) = Z(37:43);
Z_new(44) = p1_n;
Z_new(45:51) = Z(45:51);
Z_new(52) = p2_n;
Z_new(53:59) = Z(53:59);
Z_new(60) = p3_n;    
end



