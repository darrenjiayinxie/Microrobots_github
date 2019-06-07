function [F, J, domerr] = mcp_funjac_probabilistic_grasp_planning_simplified(z, jacflag)


%% initialize
z = z(:);
F = [];
J = [];
domerr = 0;

%% obtain value of global variables
global h;

global q_old ;
q_xo = q_old(1);
q_yo = q_old(2);
q_zo = q_old(3);
q0_o = q_old(4);
q1_o = q_old(5);
q2_o = q_old(6);
q3_o = q_old(7);

global nu_old;
v_xo =nu_old(1);
v_yo =nu_old(2);
v_zo =nu_old(3);
w_xo =nu_old(4);
w_yo =nu_old(5);
w_zo =nu_old(6);

global m mu_g mu_l mu_r I_xx I_yy I_zz e_o e_r e_t g;

global r heg;

global p_x p_y p_z p_xt p_yt p_zt;

% Dimension, position and velocity of gripper
global y_l_o v_l Height_l  y_r_o v_r Height_r; 

%% unknown variables
v_x = z(1);
v_y = z(2);
v_z=  z(3);
w_x = z(4);
w_y = z(5);
w_z = z(6);

a1_x_l = z(7);
a1_y_l = z(8);
a1_z_l = z(9);

a2_x_l = z(10);
a2_y_l = z(11);
a2_z_l = z(12);

a1_x_r = z(13);
a1_y_r = z(14);
a1_z_r = z(15);

a2_x_r = z(16);
a2_y_r = z(17);
a2_z_r = z(18);

a1_x_g = z(19);
a1_y_g = z(20);
a1_z_g = z(21);

a2_x_g = z(22);
a2_y_g = z(23);
a2_z_g = z(24);

p_t_l = z(25);
p_o_l = z(26);
p_r_l = z(27);

p_t_r = z(28);
p_o_r = z(29);
p_r_r = z(30);

p_t_g = z(31);
p_o_g = z(32);
p_r_g = z(33);

sig_l = z(34);
sig_r = z(35);
sig_g = z(36);

l1_l = z(37);
l2_l = z(38);
l3_l = z(39);

l1_lg = z(40);
l2_lg= z(41);

l1_r = z(42);
l2_r = z(43);
l3_r = z(44);

l1_rg = z(45);
l2_rg = z(46);

l1_g = z(47);
l2_g = z(48);
l3_g = z(49);

lg = z(50);

p_n_l = z(51);
p_n_r = z(52);
p_n_g = z(53);


%% intermediate variables appear in the chain rule
% z - vector of unknown variables 

% The intermediate variables are:
% 1. [qo(z) q1(z) q2(z) q3(z)] - Quaternions, which depends on z
% 2. Rij(qo,q1,q2,q3) - elements in rotation matrix, which depends on the quaternion
% 3. q_x(z),q_y(z),q_z(z) - position vector in l+1, which depends on q_old and z
% 4. Iij(qo,q1,q2,q3) - elements in Inerteria matrix, which depends on the quaternion

% Thus system of equations F can be described by the unknown variables z or
% with intermediate variables:
% 1. F(z,Rij,q_x,q_y,q_z) - F descriabed by z with intermediate variables
% 2. F(z) - F described by unknown variables z

% Based on the chain rule, Jocobain J = J1*J2*J3 where:
% J - dF(z)/dz Jacobian matrix 
% J1 = dF(z,Rij,Iij,q_x,q_y,q_z)/d[z,Rij,q_x,q_y,q_z] 
% J2 = d[z,Rij(qo,q1,q2,q3),q_x(z),q_y(z),q_z(z)]/d[z,qo,q1,q2,q3] 
% J3 = d[z,qo(z),q1(z),q2(z),q3(z)]/dz 

q0 = -(2*((h*q1_o*w_x)/2 - q0_o + (h*q2_o*w_y)/2 + (h*q3_o*w_z)/2))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(1/2);
q1 = (2*(q1_o + (h*q0_o*w_x)/2 + (h*q3_o*w_y)/2 - (h*q2_o*w_z)/2))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(1/2);
q2 = (2*(q2_o - (h*q3_o*w_x)/2 + (h*q0_o*w_y)/2 + (h*q1_o*w_z)/2))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(1/2);
q3 = (2*(q3_o + (h*q2_o*w_x)/2 - (h*q1_o*w_y)/2 + (h*q0_o*w_z)/2))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(1/2);

q_x = q_xo+h*v_x;
q_y = q_yo+h*v_y;
q_z = q_zo+h*v_z;

R11 = q0^2 + q1^2 - q2^2 - q3^2;
R12 = 2*q1*q2 - 2*q0*q3;
R13 = 2*q0*q2 + 2*q1*q3;
R21 = 2*q0*q3 + 2*q1*q2;
R22 = q0^2 - q1^2 + q2^2 - q3^2;
R23 = 2*q2*q3 - 2*q0*q1;
R31 = 2*q1*q3 - 2*q0*q2;
R32 = 2*q0*q1 + 2*q2*q3;
R33 = q0^2 - q1^2 - q2^2 + q3^2;

R = [R11 R12 R13; R21 R22 R23; R31 R32 R33];

r_l = [a1_x_l-q_x;a1_y_l-q_y;a1_z_l-q_z];

r_r = [a1_x_r-q_x;a1_y_r-q_y;a1_z_r-q_z];

r_g = [a1_x_g-q_x;a1_y_g-q_y;a1_z_g-q_z];

I11 = (I_xx*(q0^2 + q1^2 - q2^2 - q3^2) - I_yy*(2*q0*q3 - 2*q1*q2) + I_zz*(2*q0*q2 + 2*q1*q3))*(q0^2 + q1^2 - q2^2 - q3^2);
I12 = (2*q0*q3 + 2*q1*q2)*(I_xx*(q0^2 + q1^2 - q2^2 - q3^2) - I_yy*(2*q0*q3 - 2*q1*q2) + I_zz*(2*q0*q2 + 2*q1*q3));
I13 = -(2*q0*q2 - 2*q1*q3)*(I_xx*(q0^2 + q1^2 - q2^2 - q3^2) - I_yy*(2*q0*q3 - 2*q1*q2) + I_zz*(2*q0*q2 + 2*q1*q3));
I21 = (I_yy*(q0^2 - q1^2 + q2^2 - q3^2) + I_xx*(2*q0*q3 + 2*q1*q2) - I_zz*(2*q0*q1 - 2*q2*q3))*(q0^2 + q1^2 - q2^2 - q3^2);
I22 = (2*q0*q3 + 2*q1*q2)*(I_yy*(q0^2 - q1^2 + q2^2 - q3^2) + I_xx*(2*q0*q3 + 2*q1*q2) - I_zz*(2*q0*q1 - 2*q2*q3));
I23 = -(2*q0*q2 - 2*q1*q3)*(I_yy*(q0^2 - q1^2 + q2^2 - q3^2) + I_xx*(2*q0*q3 + 2*q1*q2) - I_zz*(2*q0*q1 - 2*q2*q3));
I31 = (I_zz*(q0^2 - q1^2 - q2^2 + q3^2) - I_xx*(2*q0*q2 - 2*q1*q3) + I_yy*(2*q0*q1 + 2*q2*q3))*(q0^2 + q1^2 - q2^2 - q3^2);
I32 = (2*q0*q3 + 2*q1*q2)*(I_zz*(q0^2 - q1^2 - q2^2 + q3^2) - I_xx*(2*q0*q2 - 2*q1*q3) + I_yy*(2*q0*q1 + 2*q2*q3));
I33 = -(2*q0*q2 - 2*q1*q3)*(I_zz*(q0^2 - q1^2 - q2^2 + q3^2) - I_xx*(2*q0*q2 - 2*q1*q3) + I_yy*(2*q0*q1 + 2*q2*q3));


%% Newton_Euler equations

F(1) = p_o_l - p_o_r + p_t_g + p_x - m*(v_x - v_xo);
F(2) = p_n_l - p_n_r + p_o_g + p_y - m*(v_y - v_yo);
F(3) = p_n_g + p_t_l + p_t_r + p_z - m*(v_z - v_zo) - g*h*m;

F(4) = p_xt - I11*(w_x - w_xo) - I12*(w_y - w_yo) - I13*(w_z - w_zo) + p_n_g*(a1_y_g - q_y) - p_n_l*(a1_z_l - q_z) + p_n_r*(a1_z_r - q_z) - p_o_g*(a1_z_g - q_z) + p_t_l*(a1_y_l - q_y) + p_t_r*(a1_y_r - q_y) - w_yo*(I31*w_xo + I32*w_yo + I33*w_zo) + w_zo*(I21*w_xo + I22*w_yo + I23*w_zo);
F(5) = p_r_l - p_r_r + p_yt - I21*(w_x - w_xo) - I22*(w_y - w_yo) - I23*(w_z - w_zo) - p_n_g*(a1_x_g - q_x) + p_o_l*(a1_z_l - q_z) - p_o_r*(a1_z_r - q_z) - p_t_l*(a1_x_l - q_x) - p_t_r*(a1_x_r - q_x) + p_t_g*(a1_z_g - q_z) + w_xo*(I31*w_xo + I32*w_yo + I33*w_zo) - w_zo*(I11*w_xo + I12*w_yo + I13*w_zo);
F(6) = p_r_g + p_zt - I31*(w_x - w_xo) - I32*(w_y - w_yo) - I33*(w_z - w_zo) + p_n_l*(a1_x_l - q_x) - p_n_r*(a1_x_r - q_x) + p_o_g*(a1_x_g - q_x) - p_o_l*(a1_y_l - q_y) + p_o_r*(a1_y_r - q_y) - p_t_g*(a1_y_g - q_y) - w_xo*(I21*w_xo + I22*w_yo + I23*w_zo) + w_yo*(I11*w_xo + I12*w_yo + I13*w_zo);

%% Friction Model without complementarity equation

F(7) = p_t_l*sig_l + e_t^2*mu_l*p_n_l*v_z - e_t^2*mu_l*p_n_l*w_y*(a1_x_l - q_x) + e_t^2*mu_l*p_n_l*w_x*(a1_y_l - q_y);
F(8) = p_o_l*sig_l + e_o^2*mu_l*p_n_l*v_x - e_o^2*mu_l*p_n_l*w_z*(a1_y_l - q_y) + e_o^2*mu_l*p_n_l*w_y*(a1_z_l - q_z);
F(9) = mu_l*p_n_l*w_y*e_r^2 + p_r_l*sig_l;

F(10) = p_t_r*sig_r + e_t^2*mu_r*p_n_r*v_z - e_t^2*mu_r*p_n_r*w_y*(a1_x_r - q_x) + e_t^2*mu_r*p_n_r*w_x*(a1_y_r - q_y);
F(11) = p_o_r*sig_r - e_o^2*mu_r*p_n_r*v_x + e_o^2*mu_r*p_n_r*w_z*(a1_y_r - q_y) - e_o^2*mu_r*p_n_r*w_y*(a1_z_r - q_z);
F(12) = - mu_r*p_n_r*w_y*e_r^2 + p_r_r*sig_r;

F(13) = p_t_g*sig_g + e_t^2*mu_g*p_n_g*v_x - e_t^2*mu_g*p_n_g*w_z*(a1_y_g - q_y) + e_t^2*mu_g*p_n_g*w_y*(a1_z_g - q_z);
F(14) = p_o_g*sig_g + e_o^2*mu_g*p_n_g*v_y + e_o^2*mu_g*p_n_g*w_z*(a1_x_g - q_x) - e_o^2*mu_g*p_n_g*w_x*(a1_z_g - q_z);
F(15) = mu_g*p_n_g*w_z*e_r^2 + p_r_g*sig_g;

%% contact constraints without complementarity equation

% left gripper
F(16) = a2_x_l - a1_x_l;
F(17) = a2_y_l - a1_y_l + l1_lg;
F(18) = a2_z_l - a1_z_l - l1_lg*l2_lg;
F(19) = R13*l3_l - R13*l2_l + l1_l*(2*R11*(R11*a1_x_l + R21*a1_y_l + R31*a1_z_l - R11*q_x - R21*q_y - R31*q_z) + 2*R12*(R12*a1_x_l + R22*a1_y_l + R32*a1_z_l - R12*q_x - R22*q_y - R32*q_z));
F(20) = R23*l3_l - R23*l2_l + l1_l*(2*R21*(R11*a1_x_l + R21*a1_y_l + R31*a1_z_l - R11*q_x - R21*q_y - R31*q_z) + 2*R22*(R12*a1_x_l + R22*a1_y_l + R32*a1_z_l - R12*q_x - R22*q_y - R32*q_z)) + 1;
F(21) = R33*l3_l - R33*l2_l - l2_lg + l1_l*(2*R31*(R11*a1_x_l + R21*a1_y_l + R31*a1_z_l - R11*q_x - R21*q_y - R31*q_z) + 2*R32*(R12*a1_x_l + R22*a1_y_l + R32*a1_z_l - R12*q_x - R22*q_y - R32*q_z));

% right gripper
F(22) = a2_x_r - a1_x_r;
F(23) = a2_y_r - a1_y_r - l1_rg;
F(24) = a2_z_r - a1_z_r - l1_rg*l2_rg;
F(25) = R13*l3_r - R13*l2_r + l1_r*(2*R11*(R11*a1_x_r + R21*a1_y_r + R31*a1_z_r - R11*q_x - R21*q_y - R31*q_z) + 2*R12*(R12*a1_x_r + R22*a1_y_r + R32*a1_z_r - R12*q_x - R22*q_y - R32*q_z));
F(26) = R23*l3_r - R23*l2_r + l1_r*(2*R21*(R11*a1_x_r + R21*a1_y_r + R31*a1_z_r - R11*q_x - R21*q_y - R31*q_z) + 2*R22*(R12*a1_x_r + R22*a1_y_r + R32*a1_z_r - R12*q_x - R22*q_y - R32*q_z)) - 1;
F(27) = R33*l3_r - R33*l2_r - l2_rg + l1_r*(2*R31*(R11*a1_x_r + R21*a1_y_r + R31*a1_z_r - R11*q_x - R21*q_y - R31*q_z) + 2*R32*(R12*a1_x_r + R22*a1_y_r + R32*a1_z_r - R12*q_x - R22*q_y - R32*q_z));

% ground 
F(28) = a2_x_g - a1_x_g;
F(29) = a2_y_g - a1_y_g;
F(30) = a2_z_g - a1_z_g + lg;
F(31) = R13*l3_g - R13*l2_g + l1_g*(2*R11*(R11*a1_x_g + R21*a1_y_g + R31*a1_z_g - R11*q_x - R21*q_y - R31*q_z) + 2*R12*(R12*a1_x_g + R22*a1_y_g + R32*a1_z_g - R12*q_x - R22*q_y - R32*q_z));
F(32) = R23*l3_g - R23*l2_g + l1_g*(2*R21*(R11*a1_x_g + R21*a1_y_g + R31*a1_z_g - R11*q_x - R21*q_y - R31*q_z) + 2*R22*(R12*a1_x_g + R22*a1_y_g + R32*a1_z_g - R12*q_x - R22*q_y - R32*q_z));
F(33) = R33*l3_g - R33*l2_g + l1_g*(2*R31*(R11*a1_x_g + R21*a1_y_g + R31*a1_z_g - R11*q_x - R21*q_y - R31*q_z) + 2*R32*(R12*a1_x_g + R22*a1_y_g + R32*a1_z_g - R12*q_x - R22*q_y - R32*q_z)) + 1;

%% complementarity equation in Friction Model
F(34) = mu_l^2*p_n_l^2 - p_r_l^2/e_r^2 - p_t_l^2/e_t^2 - p_o_l^2/e_o^2;
F(35) = mu_r^2*p_n_r^2 - p_r_r^2/e_r^2 - p_t_r^2/e_t^2 - p_o_r^2/e_o^2;
F(36) = mu_g^2*p_n_g^2 - p_r_g^2/e_r^2 - p_t_g^2/e_t^2 - p_o_g^2/e_o^2;

%% complementarity equations in contact constraints

% left gripper
F(37) = r^2 - (R12*a1_x_l + R22*a1_y_l + R32*a1_z_l - R12*q_x - R22*q_y - R32*q_z)^2 - (R11*a1_x_l + R21*a1_y_l + R31*a1_z_l - R11*q_x - R21*q_y - R31*q_z)^2;
F(38) = heg/2 + R13*a1_x_l + R23*a1_y_l + R33*a1_z_l - R13*q_x - R23*q_y - R33*q_z;
F(39) = heg/2 - R13*a1_x_l - R23*a1_y_l - R33*a1_z_l + R13*q_x + R23*q_y + R33*q_z;

F(40) =  y_l_o - a2_y_l + h*v_l;
F(41) = a2_z_l - Height_l;
        
% right gripper
F(42) = r^2 - (R12*a1_x_r + R22*a1_y_r + R32*a1_z_r - R12*q_x - R22*q_y - R32*q_z)^2 - (R11*a1_x_r + R21*a1_y_r + R31*a1_z_r - R11*q_x - R21*q_y - R31*q_z)^2;
F(43) = heg/2 + R13*a1_x_r + R23*a1_y_r + R33*a1_z_r - R13*q_x - R23*q_y - R33*q_z;
F(44) = heg/2 - R13*a1_x_r - R23*a1_y_r - R33*a1_z_r + R13*q_x + R23*q_y + R33*q_z;

F(45) = a2_y_r - y_r_o - h*v_r;
F(46) = a2_z_r - Height_r;

% ground
F(47) =  r^2 - (R12*a1_x_g + R22*a1_y_g + R32*a1_z_g - R12*q_x - R22*q_y - R32*q_z)^2 - (R11*a1_x_g + R21*a1_y_g + R31*a1_z_g - R11*q_x - R21*q_y - R31*q_z)^2;
F(48) = heg/2 + R13*a1_x_g + R23*a1_y_g + R33*a1_z_g - R13*q_x - R23*q_y - R33*q_z;
F(49) = heg/2 - R13*a1_x_g - R23*a1_y_g - R33*a1_z_g + R13*q_x + R23*q_y + R33*q_z;

F(50) = -a2_z_g;

% non-penetration
F(51) = a1_y_l - y_l_o - h*v_l;
F(52) = y_r_o - a1_y_r + h*v_r;
F(53) = a1_z_g;
%% Jacobian matrix (Using chain rules)

if (jacflag)
    
    J1 = zeros(53,74);
    
    % Newton_Euler equations
    J1(1:3,1:3) = -m*eye(3);  %(v_x v_y v_z)
    J1(1:3,[31,32,53]) = eye(3); %(p_t_g,p_o_g,p_n_g)
    J1(1:3,[26,51,25]) = eye(3); %(p_o_l,p_n_l,p_t_l)
    J1(1:3,[29,52,28]) = [-1,0,0;0,-1,0;0,0,1]; %(p_o_r,p_n_r,p_t_r)
    
    J1(4:6,4:6) = -[I11,I12,I13;I21,I22,I23;I31,I32,I33]; %(w_x w_y w_z)
    J1(4:6,7:9) = product_skew_mat([p_o_l,p_n_l,p_t_l]);  %(a1_x_l a1_y_l a1_z_l)
    J1(4:6,13:15) = product_skew_mat([p_o_r,p_n_r,-p_t_r]); %(a1_x_r a1_y_r a1_z_r)
    J1(4:6,19:21) = -product_skew_mat([p_o_g,p_n_g,p_t_g]); %(a1_x_g a1_y_g a1_z_g)
    J1(4:6,[26,51,25]) = product_skew_mat([a1_x_l - q_x,a1_y_l - q_y,a1_z_l - q_z]); %(p_o_l,p_n_l,p_t_l)
    J1(4:6,[29,52,28]) = product_skew_mat([-a1_x_r + q_x,a1_y_r - q_y,-a1_z_r + q_z]); %(p_o_r,p_n_r,p_t_r)
    J1(4:6,[31,32,53]) = product_skew_mat([a1_x_g - q_x,a1_y_g - q_y,a1_z_g - q_z]); %(p_t_g,p_o_g,p_n_g)
    J1(4:6,[27,30,33]) = [0,0,0;1,-1,0;0,0,1];  %(p_r_l,p_r_r,p_r_g)
    J1(4:6,63:71) = [ w_xo - w_x, w_yo - w_y, w_zo - w_z,  w_xo*w_zo,  w_yo*w_zo,     w_zo^2, -w_xo*w_yo,    -w_yo^2, -w_yo*w_zo   %(I11-I33)
                     -w_xo*w_zo, -w_yo*w_zo,    -w_zo^2, w_xo - w_x, w_yo - w_y, w_zo - w_z,     w_xo^2,  w_xo*w_yo,  w_xo*w_zo
                      w_xo*w_yo,     w_yo^2,  w_yo*w_zo,    -w_xo^2, -w_xo*w_yo, -w_xo*w_zo, w_xo - w_x, w_yo - w_y, w_zo - w_z];
    J1(4:6,72:74) = product_skew_mat([p_o_l - p_o_r + p_t_g, p_n_l - p_n_r + p_o_g, p_n_g + p_t_l + p_t_r]); % q_x q_y q_z
    
    % contact constraints without complementarity equation - left gripper
    J1(7:9,1:3) = [ 0, 0, e_t^2*mu_l*p_n_l;e_o^2*mu_l*p_n_l, 0, 0; 0, 0, 0];     %v_x v_y v_z
    J1(7:9,4:6) = [ e_t^2*mu_l*p_n_l*(a1_y_l - q_y), -e_t^2*mu_l*p_n_l*(a1_x_l - q_x),                                0  % w_x w_y w_z
                                                  0,  e_o^2*mu_l*p_n_l*(a1_z_l - q_z), -e_o^2*mu_l*p_n_l*(a1_y_l - q_y)
                                                  0,                 e_r^2*mu_l*p_n_l,                                0]; 
    J1(7:9,7:9) = [ -e_t^2*mu_l*p_n_l*w_y,  e_t^2*mu_l*p_n_l*w_x,                    0    %a1_x_l a1_y_l a1_z_l
                                        0, -e_o^2*mu_l*p_n_l*w_z, e_o^2*mu_l*p_n_l*w_y
                                        0,                     0,                    0];  
    J1(7:9,25:27) = sig_l*eye(3);  %p_t_l p_o_l p_r_l
    J1(7:9,34) = [p_t_l;p_o_l;p_r_l]; % sig_l
    J1(7:9,51) = [e_t^2*mu_l*v_z - e_t^2*mu_l*w_y*(a1_x_l - q_x) + e_t^2*mu_l*w_x*(a1_y_l - q_y);e_o^2*mu_l*v_x - e_o^2*mu_l*w_z*(a1_y_l - q_y) + e_o^2*mu_l*w_y*(a1_z_l - q_z);e_r^2*mu_l*w_y]; %p_n_l
    J1(7:9,72:74) = [ e_t^2*mu_l*p_n_l*w_y, -e_t^2*mu_l*p_n_l*w_x,                     0
                                         0,  e_o^2*mu_l*p_n_l*w_z, -e_o^2*mu_l*p_n_l*w_y
                                         0,                     0,                     0];  % q_x q_y q_z
    % contact constraints without complementarity equation - right gripper
    J1(10:12,1:3) = [ 0, 0, e_t^2*mu_r*p_n_r;-e_o^2*mu_r*p_n_r, 0, 0; 0, 0, 0];     %v_x v_y v_z
    J1(10:12,4:6) = [ e_t^2*mu_r*p_n_r*(a1_y_r - q_y), -e_t^2*mu_r*p_n_r*(a1_x_r - q_x),                                0  % w_x w_y w_z
                                                  0,   -e_o^2*mu_r*p_n_r*(a1_z_r - q_z),  e_o^2*mu_r*p_n_r*(a1_y_r - q_y)
                                                  0,                  -e_r^2*mu_r*p_n_r,                                0]; 
    J1(10:12,13:15) = [ -e_t^2*mu_r*p_n_r*w_y,  e_t^2*mu_r*p_n_r*w_x,                    0    %a1_x_r a1_y_r a1_z_r
                                        0,    e_o^2*mu_r*p_n_r*w_z, -e_o^2*mu_r*p_n_r*w_y
                                        0,                     0,                    0];  
    J1(10:12,28:30) = sig_r*eye(3);  %p_t_r p_o_r p_r_r
    J1(10:12,35) = [p_t_r;p_o_r;p_r_r]; % sig_r
    J1(10:12,52) = [e_t^2*mu_r*v_z - e_t^2*mu_r*w_y*(a1_x_r - q_x) + e_t^2*mu_r*w_x*(a1_y_r - q_y);
                    e_o^2*mu_r*w_z*(a1_y_r - q_y) - e_o^2*mu_r*v_x - e_o^2*mu_r*w_y*(a1_z_r - q_z);
                    -e_r^2*mu_r*w_y]; %p_n_r
    J1(10:12,72:74) = [ e_t^2*mu_r*p_n_r*w_y, -e_t^2*mu_r*p_n_r*w_x,                    0
                                           0, -e_o^2*mu_r*p_n_r*w_z, e_o^2*mu_r*p_n_r*w_y
                                           0,                     0,                    0];  % q_x q_y q_z
                                       
    % contact constraints without complementarity equation - ground
    J1(13:15,1:3) = [ e_t^2*mu_g*p_n_g, 0, 0; 0, e_o^2*mu_g*p_n_g, 0; 0, 0, 0];  %v_x v_y v_z
    J1(13:15,4:6) = [ 0, e_t^2*mu_g*p_n_g*(a1_z_g - q_z), -e_t^2*mu_g*p_n_g*(a1_y_g - q_y) %w_x w_y w_z
                     -e_o^2*mu_g*p_n_g*(a1_z_g - q_z), 0,  e_o^2*mu_g*p_n_g*(a1_x_g - q_x)
                                                    0, 0,                e_r^2*mu_g*p_n_g];
    J1(13:15,19:21) = [                    0, -e_t^2*mu_g*p_n_g*w_z,  e_t^2*mu_g*p_n_g*w_y
                       e_o^2*mu_g*p_n_g*w_z,                     0, -e_o^2*mu_g*p_n_g*w_x
                                          0,                     0,                     0];  %a1_x_g a1_y_g a1_z_g
    J1(13:15,31:33) = sig_g*eye(3); % p_t_g p_o_g p_r_g
    J1(13:15,36) = [p_t_g;p_o_g;p_r_g]; % sig_g
    J1(13:15,53) = [e_t^2*mu_g*v_x - e_t^2*mu_g*w_z*(a1_y_g - q_y) + e_t^2*mu_g*w_y*(a1_z_g - q_z)
                    e_o^2*mu_g*v_y + e_o^2*mu_g*w_z*(a1_x_g - q_x) - e_o^2*mu_g*w_x*(a1_z_g - q_z)
                                                                                   e_r^2*mu_g*w_z]; %p_n_g
    J1(13:15,72:74) = [                     0, e_t^2*mu_g*p_n_g*w_z, -e_t^2*mu_g*p_n_g*w_y
                        -e_o^2*mu_g*p_n_g*w_z,                    0,  e_o^2*mu_g*p_n_g*w_x
                                            0,                    0,                     0]; % q_x q_y q_z
                                                                               
    % contact constraints without complementarity equation - left gripper
    J1(16:18,7:9) = -eye(3); % a1_x_l a1_y_l a1_z_l
    J1(16:18,10:12) = eye(3); % a2_x_l a2_y_l a2_z_l
    J1(16:18,40) = [0;1;-l2_lg]; % l1_lg
    J1(16:18,41) = [0;0;-l1_lg]; % l2_lg
    J1(19:21,7:9) = 2*l1_l*R*(R'-[0,0,0;0,0,0;R13,R23,R33]); % a1_x_l a1_y_l a1_z_l
    J1(19:21,37:39) = [2*(R-[0,0,R13;0,0,R23;0,0,R33])*R'*r_l,-R(:,3),R(:,3)]; % l1_l l2_l l3_l
    J1(19:21,41) = [0;0;-1]; %l2_lg
    J1(19:21,[54,57,60]) = 2*l1_l*(R(:,1)*r_l'+R(:,1)'*r_l*eye(3)); %R11 R21 R31
    J1(19:21,[55,58,61]) = 2*l1_l*(R(:,2)*r_l'+R(:,2)'*r_l*eye(3)); %R12 R22 R32
    J1(19:21,[56,59,62]) = (l3_l - l2_l)*eye(3);                    %R13 R23 R33
    J1(19:21,72:74) = -2*l1_l*R*(R'-[0,0,0;0,0,0;R13,R23,R33]);     %q_x q_y q_z
    
    % contact constraints without complementarity equation - right gripper
    J1(22:24,13:15) = -eye(3);% a1_x_r a1_y_r a1_z_r
    J1(22:24,16:18) = eye(3); % a2_x_r a2_y_r a2_z_r
    J1(22:24,45) = [0;-1;-l2_rg]; % l1_rg
    J1(22:24,46) = [0;0;-l1_rg]; % l2_rg
    J1(25:27,13:15) = 2*l1_r*R*(R'-[0,0,0;0,0,0;R13,R23,R33]); % a1_x_r a1_y_r a1_z_r
    J1(25:27,42:44) = [2*(R-[0,0,R13;0,0,R23;0,0,R33])*R'*r_r,-R(:,3),R(:,3)]; % l1_r l2_r l3_r
    J1(25:27,46) = [0;0;-1]; %l2_rg
    J1(25:27,[54,57,60]) = 2*l1_r*(R(:,1)*r_r'+R(:,1)'*r_r*eye(3)); %R11 R21 R31
    J1(25:27,[55,58,61]) = 2*l1_r*(R(:,2)*r_r'+R(:,2)'*r_r*eye(3)); %R12 R22 R32
    J1(25:27,[56,59,62]) = (l3_r - l2_r)*eye(3);                    %R13 R23 R33
    J1(25:27,72:74) = -2*l1_r*R*(R'-[0,0,0;0,0,0;R13,R23,R33]);     %q_x q_y q_z
    
    % contact constraints without complementarity equation - ground
    J1(28:30,19:21) = -eye(3); % a1_x_g a1_y_g a1_z_g 
    J1(28:30,22:24) =  eye(3); % a2_x_g a2_y_g a2_z_g
    J1(28:30,50) = [0;0;1]; % l_g
    J1(31:33,19:21) = 2*l1_g*R*(R'-[0,0,0;0,0,0;R13,R23,R33]); % a1_x_g a1_y_g a1_z_g
    J1(31:33,47:49) = [2*(R-[0,0,R13;0,0,R23;0,0,R33])*R'*r_g,-R(:,3),R(:,3)]; % l1_g l2_g l3_g
    J1(31:33,[54,57,60]) =  2*l1_g*(R(:,1)*r_g'+R(:,1)'*r_g*eye(3)); %R11 R21 R31
    J1(31:33,[55,58,61]) = 2*l1_g*(R(:,2)*r_g'+R(:,2)'*r_g*eye(3)); %R12 R22 R32
    J1(31:33,[56,59,62]) = (l3_g - l2_g)*eye(3);                    %R13 R23 R33
    J1(31:33,72:74) = -2*l1_g*R*(R'-[0,0,0;0,0,0;R13,R23,R33]);     %q_x q_y q_z
    
    
    % complementarity equation in Friction Model
    J1(34,25:27) = [-(2*p_t_l)/e_t^2, -(2*p_o_l)/e_o^2, -(2*p_r_l)/e_r^2];
    J1(35,28:30) = [-(2*p_t_r)/e_t^2, -(2*p_o_r)/e_o^2, -(2*p_r_r)/e_r^2];
    J1(36,31:33) = [-(2*p_t_g)/e_t^2, -(2*p_o_g)/e_o^2, -(2*p_r_g)/e_r^2];
    J1(34:36,51:53) = [ 2*mu_l^2*p_n_l,              0,              0
                                     0, 2*mu_r^2*p_n_r,              0
                                     0,              0, 2*mu_g^2*p_n_g];
                                 
    % complementarity equations in contact constraints - left gripper
    J1(37:39,7:9) = [-2*(R-[0,0,R13;0,0,R23;0,0,R33])*R'*r_l,R(:,3),-R(:,3)]'; % a1_x_l a1_y_l a1_z_l
    temp_l = -2*R'*r_l*r_l';
    J1(37,[54,57,60]) = temp_l(1,:);    %R11 R21 R31
    J1(37,[55,58,61]) = temp_l(2,:);    %R12 R22 R32
    J1(37:39,[56,59,62]) = [0,0,0;r_l';-r_l'];   %R13 R23 R33
    J1(37:39,72:74) = [2*(R-[0,0,R13;0,0,R23;0,0,R33])*R'*r_l,-R(:,3),R(:,3)]';   %q_x q_y q_z
    J1(40:41,10:12) = [ 0, -1, 0 
                        0,  0, 1];  % a2_x_l a2_y_l a2_z_l
                         
     % complementarity equations in contact constraints - right gripper
    J1(42:44,13:15) = [-2*(R-[0,0,R13;0,0,R23;0,0,R33])*R'*r_r,R(:,3),-R(:,3)]'; % a1_x_r a1_y_r a1_z_r
    temp_r = -2*R'*r_r*r_r';
    J1(42,[54,57,60]) = temp_r(1,:);    %R11 R21 R31
    J1(42,[55,58,61]) = temp_r(2,:);    %R12 R22 R32
    J1(42:44,[56,59,62]) = [0,0,0;r_r';-r_r'];   %R13 R23 R33
    J1(42:44,72:74) = [2*(R-[0,0,R13;0,0,R23;0,0,R33])*R'*r_r,-R(:,3),R(:,3)]';   %q_x q_y q_z
    J1(45:46,16:18) = [ 0, 1, 0
                        0, 0, 1];  % a2_x_r a2_y_r a2_z_r
    
                    
     % complementarity equations in contact constraints - ground
    J1(47:49,19:21) = [-2*(R-[0,0,R13;0,0,R23;0,0,R33])*R'*r_g,R(:,3),-R(:,3)]'; % a1_x_g a1_y_g a1_z_g
    temp_g = -2*R'*r_g*r_g';
    J1(47,[54,57,60]) = temp_g(1,:);    %R11 R21 R31
    J1(47,[55,58,61]) = temp_g(2,:);    %R12 R22 R32
    J1(47:49,[56,59,62]) = [0,0,0;r_g';-r_g'];   %R13 R23 R33
    J1(47:49,72:74) = [2*(R-[0,0,R13;0,0,R23;0,0,R33])*R'*r_g,-R(:,3),R(:,3)]';   %q_x q_y q_z
    J1(50,24) = -1;                
     % non-penetration
    J1(51:53,[8,14,21]) = [1,0,0;0,-1,0;0,0,1];
    
    J2 = zeros(74,57);
    J2(1:53,1:53) = eye(53);
    J2(54:62,end-3:end) = [2*q0,  2*q1, -2*q2, -2*q3
                          -2*q3,  2*q2,  2*q1, -2*q0
                           2*q2,  2*q3,  2*q0,  2*q1
                           2*q3,  2*q2,  2*q1,  2*q0
                           2*q0, -2*q1,  2*q2, -2*q3
                          -2*q1, -2*q0,  2*q3,  2*q2
                          -2*q2,  2*q3, -2*q0,  2*q1
                           2*q1,  2*q0,  2*q3,  2*q2
                           2*q0, -2*q1, -2*q2,  2*q3];
    J2(63:71,end-3:end) = [ 2*q0*(I_xx*(q0^2 + q1^2 - q2^2 - q3^2) - I_yy*(2*q0*q3 - 2*q1*q2) + I_zz*(2*q0*q2 + 2*q1*q3)) + (2*I_xx*q0 - 2*I_yy*q3 + 2*I_zz*q2)*(q0^2 + q1^2 - q2^2 - q3^2), 2*q1*(I_xx*(q0^2 + q1^2 - q2^2 - q3^2) - I_yy*(2*q0*q3 - 2*q1*q2) + I_zz*(2*q0*q2 + 2*q1*q3)) + (2*I_xx*q1 + 2*I_yy*q2 + 2*I_zz*q3)*(q0^2 + q1^2 - q2^2 - q3^2),   (2*I_yy*q1 - 2*I_xx*q2 + 2*I_zz*q0)*(q0^2 + q1^2 - q2^2 - q3^2) - 2*q2*(I_xx*(q0^2 + q1^2 - q2^2 - q3^2) - I_yy*(2*q0*q3 - 2*q1*q2) + I_zz*(2*q0*q2 + 2*q1*q3)), - 2*q3*(I_xx*(q0^2 + q1^2 - q2^2 - q3^2) - I_yy*(2*q0*q3 - 2*q1*q2) + I_zz*(2*q0*q2 + 2*q1*q3)) - (2*I_xx*q3 + 2*I_yy*q0 - 2*I_zz*q1)*(q0^2 + q1^2 - q2^2 - q3^2)
                                   2*q3*(I_xx*(q0^2 + q1^2 - q2^2 - q3^2) - I_yy*(2*q0*q3 - 2*q1*q2) + I_zz*(2*q0*q2 + 2*q1*q3)) + (2*q0*q3 + 2*q1*q2)*(2*I_xx*q0 - 2*I_yy*q3 + 2*I_zz*q2),         2*q2*(I_xx*(q0^2 + q1^2 - q2^2 - q3^2) - I_yy*(2*q0*q3 - 2*q1*q2) + I_zz*(2*q0*q2 + 2*q1*q3)) + (2*q0*q3 + 2*q1*q2)*(2*I_xx*q1 + 2*I_yy*q2 + 2*I_zz*q3),           2*q1*(I_xx*(q0^2 + q1^2 - q2^2 - q3^2) - I_yy*(2*q0*q3 - 2*q1*q2) + I_zz*(2*q0*q2 + 2*q1*q3)) + (2*q0*q3 + 2*q1*q2)*(2*I_yy*q1 - 2*I_xx*q2 + 2*I_zz*q0),           2*q0*(I_xx*(q0^2 + q1^2 - q2^2 - q3^2) - I_yy*(2*q0*q3 - 2*q1*q2) + I_zz*(2*q0*q2 + 2*q1*q3)) - (2*q0*q3 + 2*q1*q2)*(2*I_xx*q3 + 2*I_yy*q0 - 2*I_zz*q1)
                                 - 2*q2*(I_xx*(q0^2 + q1^2 - q2^2 - q3^2) - I_yy*(2*q0*q3 - 2*q1*q2) + I_zz*(2*q0*q2 + 2*q1*q3)) - (2*q0*q2 - 2*q1*q3)*(2*I_xx*q0 - 2*I_yy*q3 + 2*I_zz*q2),         2*q3*(I_xx*(q0^2 + q1^2 - q2^2 - q3^2) - I_yy*(2*q0*q3 - 2*q1*q2) + I_zz*(2*q0*q2 + 2*q1*q3)) - (2*q0*q2 - 2*q1*q3)*(2*I_xx*q1 + 2*I_yy*q2 + 2*I_zz*q3),         - 2*q0*(I_xx*(q0^2 + q1^2 - q2^2 - q3^2) - I_yy*(2*q0*q3 - 2*q1*q2) + I_zz*(2*q0*q2 + 2*q1*q3)) - (2*q0*q2 - 2*q1*q3)*(2*I_yy*q1 - 2*I_xx*q2 + 2*I_zz*q0),           2*q1*(I_xx*(q0^2 + q1^2 - q2^2 - q3^2) - I_yy*(2*q0*q3 - 2*q1*q2) + I_zz*(2*q0*q2 + 2*q1*q3)) + (2*q0*q2 - 2*q1*q3)*(2*I_xx*q3 + 2*I_yy*q0 - 2*I_zz*q1)
                            2*q0*(I_yy*(q0^2 - q1^2 + q2^2 - q3^2) + I_xx*(2*q0*q3 + 2*q1*q2) - I_zz*(2*q0*q1 - 2*q2*q3)) + (2*I_xx*q3 + 2*I_yy*q0 - 2*I_zz*q1)*(q0^2 + q1^2 - q2^2 - q3^2), 2*q1*(I_yy*(q0^2 - q1^2 + q2^2 - q3^2) + I_xx*(2*q0*q3 + 2*q1*q2) - I_zz*(2*q0*q1 - 2*q2*q3)) - (2*I_yy*q1 - 2*I_xx*q2 + 2*I_zz*q0)*(q0^2 + q1^2 - q2^2 - q3^2),   (2*I_xx*q1 + 2*I_yy*q2 + 2*I_zz*q3)*(q0^2 + q1^2 - q2^2 - q3^2) - 2*q2*(I_yy*(q0^2 - q1^2 + q2^2 - q3^2) + I_xx*(2*q0*q3 + 2*q1*q2) - I_zz*(2*q0*q1 - 2*q2*q3)),   (2*I_xx*q0 - 2*I_yy*q3 + 2*I_zz*q2)*(q0^2 + q1^2 - q2^2 - q3^2) - 2*q3*(I_yy*(q0^2 - q1^2 + q2^2 - q3^2) + I_xx*(2*q0*q3 + 2*q1*q2) - I_zz*(2*q0*q1 - 2*q2*q3))
                                   2*q3*(I_yy*(q0^2 - q1^2 + q2^2 - q3^2) + I_xx*(2*q0*q3 + 2*q1*q2) - I_zz*(2*q0*q1 - 2*q2*q3)) + (2*q0*q3 + 2*q1*q2)*(2*I_xx*q3 + 2*I_yy*q0 - 2*I_zz*q1),         2*q2*(I_yy*(q0^2 - q1^2 + q2^2 - q3^2) + I_xx*(2*q0*q3 + 2*q1*q2) - I_zz*(2*q0*q1 - 2*q2*q3)) - (2*q0*q3 + 2*q1*q2)*(2*I_yy*q1 - 2*I_xx*q2 + 2*I_zz*q0),           2*q1*(I_yy*(q0^2 - q1^2 + q2^2 - q3^2) + I_xx*(2*q0*q3 + 2*q1*q2) - I_zz*(2*q0*q1 - 2*q2*q3)) + (2*q0*q3 + 2*q1*q2)*(2*I_xx*q1 + 2*I_yy*q2 + 2*I_zz*q3),           2*q0*(I_yy*(q0^2 - q1^2 + q2^2 - q3^2) + I_xx*(2*q0*q3 + 2*q1*q2) - I_zz*(2*q0*q1 - 2*q2*q3)) + (2*q0*q3 + 2*q1*q2)*(2*I_xx*q0 - 2*I_yy*q3 + 2*I_zz*q2)
                                 - 2*q2*(I_yy*(q0^2 - q1^2 + q2^2 - q3^2) + I_xx*(2*q0*q3 + 2*q1*q2) - I_zz*(2*q0*q1 - 2*q2*q3)) - (2*q0*q2 - 2*q1*q3)*(2*I_xx*q3 + 2*I_yy*q0 - 2*I_zz*q1),         2*q3*(I_yy*(q0^2 - q1^2 + q2^2 - q3^2) + I_xx*(2*q0*q3 + 2*q1*q2) - I_zz*(2*q0*q1 - 2*q2*q3)) + (2*q0*q2 - 2*q1*q3)*(2*I_yy*q1 - 2*I_xx*q2 + 2*I_zz*q0),         - 2*q0*(I_yy*(q0^2 - q1^2 + q2^2 - q3^2) + I_xx*(2*q0*q3 + 2*q1*q2) - I_zz*(2*q0*q1 - 2*q2*q3)) - (2*q0*q2 - 2*q1*q3)*(2*I_xx*q1 + 2*I_yy*q2 + 2*I_zz*q3),           2*q1*(I_yy*(q0^2 - q1^2 + q2^2 - q3^2) + I_xx*(2*q0*q3 + 2*q1*q2) - I_zz*(2*q0*q1 - 2*q2*q3)) - (2*q0*q2 - 2*q1*q3)*(2*I_xx*q0 - 2*I_yy*q3 + 2*I_zz*q2)
                            2*q0*(I_zz*(q0^2 - q1^2 - q2^2 + q3^2) - I_xx*(2*q0*q2 - 2*q1*q3) + I_yy*(2*q0*q1 + 2*q2*q3)) + (2*I_yy*q1 - 2*I_xx*q2 + 2*I_zz*q0)*(q0^2 + q1^2 - q2^2 - q3^2), 2*q1*(I_zz*(q0^2 - q1^2 - q2^2 + q3^2) - I_xx*(2*q0*q2 - 2*q1*q3) + I_yy*(2*q0*q1 + 2*q2*q3)) + (2*I_xx*q3 + 2*I_yy*q0 - 2*I_zz*q1)*(q0^2 + q1^2 - q2^2 - q3^2), - 2*q2*(I_zz*(q0^2 - q1^2 - q2^2 + q3^2) - I_xx*(2*q0*q2 - 2*q1*q3) + I_yy*(2*q0*q1 + 2*q2*q3)) - (2*I_xx*q0 - 2*I_yy*q3 + 2*I_zz*q2)*(q0^2 + q1^2 - q2^2 - q3^2),   (2*I_xx*q1 + 2*I_yy*q2 + 2*I_zz*q3)*(q0^2 + q1^2 - q2^2 - q3^2) - 2*q3*(I_zz*(q0^2 - q1^2 - q2^2 + q3^2) - I_xx*(2*q0*q2 - 2*q1*q3) + I_yy*(2*q0*q1 + 2*q2*q3))
                                   2*q3*(I_zz*(q0^2 - q1^2 - q2^2 + q3^2) - I_xx*(2*q0*q2 - 2*q1*q3) + I_yy*(2*q0*q1 + 2*q2*q3)) + (2*q0*q3 + 2*q1*q2)*(2*I_yy*q1 - 2*I_xx*q2 + 2*I_zz*q0),         2*q2*(I_zz*(q0^2 - q1^2 - q2^2 + q3^2) - I_xx*(2*q0*q2 - 2*q1*q3) + I_yy*(2*q0*q1 + 2*q2*q3)) + (2*q0*q3 + 2*q1*q2)*(2*I_xx*q3 + 2*I_yy*q0 - 2*I_zz*q1),           2*q1*(I_zz*(q0^2 - q1^2 - q2^2 + q3^2) - I_xx*(2*q0*q2 - 2*q1*q3) + I_yy*(2*q0*q1 + 2*q2*q3)) - (2*q0*q3 + 2*q1*q2)*(2*I_xx*q0 - 2*I_yy*q3 + 2*I_zz*q2),           2*q0*(I_zz*(q0^2 - q1^2 - q2^2 + q3^2) - I_xx*(2*q0*q2 - 2*q1*q3) + I_yy*(2*q0*q1 + 2*q2*q3)) + (2*q0*q3 + 2*q1*q2)*(2*I_xx*q1 + 2*I_yy*q2 + 2*I_zz*q3)
                                 - 2*q2*(I_zz*(q0^2 - q1^2 - q2^2 + q3^2) - I_xx*(2*q0*q2 - 2*q1*q3) + I_yy*(2*q0*q1 + 2*q2*q3)) - (2*q0*q2 - 2*q1*q3)*(2*I_yy*q1 - 2*I_xx*q2 + 2*I_zz*q0),         2*q3*(I_zz*(q0^2 - q1^2 - q2^2 + q3^2) - I_xx*(2*q0*q2 - 2*q1*q3) + I_yy*(2*q0*q1 + 2*q2*q3)) - (2*q0*q2 - 2*q1*q3)*(2*I_xx*q3 + 2*I_yy*q0 - 2*I_zz*q1),           (2*q0*q2 - 2*q1*q3)*(2*I_xx*q0 - 2*I_yy*q3 + 2*I_zz*q2) - 2*q0*(I_zz*(q0^2 - q1^2 - q2^2 + q3^2) - I_xx*(2*q0*q2 - 2*q1*q3) + I_yy*(2*q0*q1 + 2*q2*q3)),           2*q1*(I_zz*(q0^2 - q1^2 - q2^2 + q3^2) - I_xx*(2*q0*q2 - 2*q1*q3) + I_yy*(2*q0*q1 + 2*q2*q3)) - (2*q0*q2 - 2*q1*q3)*(2*I_xx*q1 + 2*I_yy*q2 + 2*I_zz*q3)];
    J2(72:74,1:3) = h*eye(3);                   
    
    J3 = zeros(57,53);
    J3(1:53,1:53) = eye(53);
    J3(54:57,4:6) = [ -(h*(q1_o*h^2*w_y^2 - q2_o*w_x*h^2*w_y + q1_o*h^2*w_z^2 - q3_o*w_x*h^2*w_z + 2*q0_o*w_x*h + 4*q1_o))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(3/2), -(h*(q2_o*h^2*w_x^2 - q1_o*w_y*h^2*w_x + q2_o*h^2*w_z^2 - q3_o*w_y*h^2*w_z + 2*q0_o*w_y*h + 4*q2_o))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(3/2), -(h*(q3_o*h^2*w_x^2 - q1_o*w_z*h^2*w_x + q3_o*h^2*w_y^2 - q2_o*w_z*h^2*w_y + 2*q0_o*w_z*h + 4*q3_o))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(3/2)
                       (h*(q0_o*h^2*w_y^2 - q3_o*w_x*h^2*w_y + q0_o*h^2*w_z^2 + q2_o*w_x*h^2*w_z - 2*q1_o*w_x*h + 4*q0_o))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(3/2),  (h*(q3_o*h^2*w_x^2 - q0_o*w_y*h^2*w_x + q3_o*h^2*w_z^2 + q2_o*w_y*h^2*w_z - 2*q1_o*w_y*h + 4*q3_o))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(3/2), -(h*(q2_o*h^2*w_x^2 + q0_o*w_z*h^2*w_x + q2_o*h^2*w_y^2 + q3_o*w_z*h^2*w_y + 2*q1_o*w_z*h + 4*q2_o))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(3/2)
                      -(h*(q3_o*h^2*w_y^2 + q0_o*w_x*h^2*w_y + q3_o*h^2*w_z^2 + q1_o*w_x*h^2*w_z + 2*q2_o*w_x*h + 4*q3_o))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(3/2),  (h*(q0_o*h^2*w_x^2 + q3_o*w_y*h^2*w_x + q0_o*h^2*w_z^2 - q1_o*w_y*h^2*w_z - 2*q2_o*w_y*h + 4*q0_o))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(3/2),  (h*(q1_o*h^2*w_x^2 + q3_o*w_z*h^2*w_x + q1_o*h^2*w_y^2 - q0_o*w_z*h^2*w_y - 2*q2_o*w_z*h + 4*q1_o))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(3/2)
                       (h*(q2_o*h^2*w_y^2 + q1_o*w_x*h^2*w_y + q2_o*h^2*w_z^2 - q0_o*w_x*h^2*w_z - 2*q3_o*w_x*h + 4*q2_o))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(3/2), -(h*(q1_o*h^2*w_x^2 + q2_o*w_y*h^2*w_x + q1_o*h^2*w_z^2 + q0_o*w_y*h^2*w_z + 2*q3_o*w_y*h + 4*q1_o))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(3/2),  (h*(q0_o*h^2*w_x^2 - q2_o*w_z*h^2*w_x + q0_o*h^2*w_y^2 + q1_o*w_z*h^2*w_y - 2*q3_o*w_z*h + 4*q0_o))/(h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(3/2)];
 
    J = J1*J2*J3;
    J = sparse(J);
end
end
