function [F, J, domerr] = mcp_funjac_microrobot_spiked_end_tipping(z, jacflag)

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

global m mu I_xx I_yy I_zz e_o e_r e_t g;

global len1 len2 wid1 wid2 heg phi;

global p_x p_y p_z p_xt p_yt p_zt;

global M_r B F_elect V_m;

global fqn time;

global Van theta;
%% unknown variables

v_x = z(1);
v_y = z(2);
v_z = z(3);
w_x = z(4);
w_y = z(5);
w_z = z(6);

a1_x = z(7);
a1_y = z(8);
a1_z = z(9);

a2_x = z(10);
a2_y = z(11);
a2_z = z(12);

p_t = z(13);
p_o = z(14);
p_r = z(15);

sig = z(16);

l1 = z(17);
l2 = z(18);
l3 = z(19);
l4 = z(20);
l5 = z(21);
l6 = z(22);
l7 = z(23);
l8 = z(24);
l9 = z(25);
l10 = z(26);
l11 = z(27);

p_n = z(28);


%% intermediate variables appear in the chain rule
% z - vector of unknown variables 

% The intermediate variables are:
% 1. [qo(z) q1(z) q2(z) q3(z)] - Quaternions, which depends on z
% 2. Rij(qo,q1,q2,q3) - elements in rotation matrix, which depends on the quaternion
% 3. Iij(qo,q1,q2,q3) - elements in RIR', which depends on the quaternion
% 3. q_x(z),q_y(z),q_z(z) - position vector in l+1, which depends on q_old and z

% Thus system of equations F can be described by the unknown variables z or
% with intermediate variables:
% 1. F(z,Rij,Iij,q_x,q_y,q_z) - F descriabed by z with intermediate variables
% 2. F(z) - F described by unknown variables z

% Based on the chain rule, Jocobain J = J1*J2*J3 where:
% J - dF(z)/dz Jacobian matrix 
% J1 = dF(z,Rij,Iij,q_x,q_y,q_z)/d[z,Rij,q_x,q_y,q_z] 
% J2 = d[z,Rij(qo,q1,q2,q3),Iij(qo,q1,q2,q3),q_x(z),q_y(z),q_z(z)]/d[z,qo,q1,q2,q3] 
% J3 = d[z,qo(z),q1(z),q2(z),q3(z)]/dz 

%% intermediate variables
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

I11 = I_xx*R11^2 + I_yy*R12^2 + I_zz*R13^2;
I12 = I_xx*R21*R11 + I_yy* R22*R12 + I_zz*R23*R13;
I13 = I_xx*R31*R11 + I_yy* R32*R12 + I_zz*R33*R13;
I21 = I12;
I22 = I_xx*R21^2 + I_yy*R22^2 + I_zz*R23^2;
I23 = I_xx*R31*R21 + I_yy*R32*R22 + I_zz*R33*R23;
I31 = I13;
I32 = I23;
I33 = I_xx*R31^2 + I_yy*R32^2 + I_zz*R33^2;

%% Newton_Euler equations
F(1) = p_t + p_x - m*(v_x - v_xo) - g*h*m*sin(theta);
F(2) = p_o + p_y - m*(v_y - v_yo);
F(3) = p_n + p_z - F_elect*h - m*(v_z - v_zo) - g*h*m*cos(theta);

F(4) = p_xt - h*(w_yo*(I31*w_xo + I32*w_yo + I33*w_zo) - w_zo*(I21*w_xo + I22*w_yo + I23*w_zo)) - I11*(w_x - w_xo) - I12*(w_y - w_yo) - I13*(w_z - w_zo) + p_n*(a1_y - q_y) - p_o*(a1_z - q_z) + V_m*h*(B*M_r*cos(2*pi*fqn*time)*(R32*cos(phi) + R33*sin(phi)) - B*M_r*sin(2*pi*fqn*time)*(R22*cos(phi) + R23*sin(phi)));
F(5) = p_yt + h*(w_xo*(I31*w_xo + I32*w_yo + I33*w_zo) - w_zo*(I11*w_xo + I12*w_yo + I13*w_zo)) - I21*(w_x - w_xo) - I22*(w_y - w_yo) - I23*(w_z - w_zo) - p_n*(a1_x - q_x) + p_t*(a1_z - q_z) + B*M_r*V_m*h*sin(2*pi*fqn*time)*(R12*cos(phi) + R13*sin(phi));
F(6) = p_r + p_zt - h*(w_xo*(I21*w_xo + I22*w_yo + I23*w_zo) - w_yo*(I11*w_xo + I12*w_yo + I13*w_zo)) - I31*(w_x - w_xo) - I32*(w_y - w_yo) - I33*(w_z - w_zo) + p_o*(a1_x - q_x) - p_t*(a1_y - q_y) - B*M_r*V_m*h*cos(2*pi*fqn*time)*(R12*cos(phi) + R13*sin(phi));

%% Friction Model without complementarity equation

F(7) = p_t*sig + e_t^2*mu*v_x*(Van + p_n) - e_t^2*mu*w_z*(Van + p_n)*(a1_y - q_y) + e_t^2*mu*w_y*(Van + p_n)*(a1_z - q_z);
F(8) = p_o*sig + e_o^2*mu*v_y*(Van + p_n) + e_o^2*mu*w_z*(Van + p_n)*(a1_x - q_x) - e_o^2*mu*w_x*(Van + p_n)*(a1_z - q_z);
F(9) = mu*w_z*(Van + p_n)*e_r^2 + p_r*sig;

%% contact constraints without complementarity equation

F(10) = a2_x - a1_x;
F(11) = a2_y - a1_y;
F(12) = a2_z - a1_z + l11;

F(13) = R11*l1 - R11*l2 + R12*l3 - R12*l4 + R13*l9 - R13*l10 + l5*(R13 + (2*R12*wid2)/len2) - l6*(R13 + (2*R12*wid2)/len2) - l7*(R13 - (2*R12*wid2)/len2) + l8*(R13 - (2*R12*wid2)/len2);
F(14) = R21*l1 - R21*l2 + R22*l3 - R22*l4 + R23*l9 - R23*l10 + l5*(R23 + (2*R22*wid2)/len2) - l6*(R23 + (2*R22*wid2)/len2) - l7*(R23 - (2*R22*wid2)/len2) + l8*(R23 - (2*R22*wid2)/len2);
F(15) = R31*l1 - R31*l2 + R32*l3 - R32*l4 + R33*l9 - R33*l10 + l5*(R33 + (2*R32*wid2)/len2) - l6*(R33 + (2*R32*wid2)/len2) - l7*(R33 - (2*R32*wid2)/len2) + l8*(R33 - (2*R32*wid2)/len2) + 1;

%% complementarity equation in Friction Model
F(16) = mu^2*(Van + p_n)^2 - p_o^2/e_o^2 - p_r^2/e_r^2 - p_t^2/e_t^2;

%% complementarity equations in contact constraints
F(17) = heg/2 - R11*a1_x - R21*a1_y - R31*a1_z + R11*q_x + R21*q_y + R31*q_z;
F(18) = heg/2 + R11*a1_x + R21*a1_y + R31*a1_z - R11*q_x - R21*q_y - R31*q_z;
F(19) = len1/2 - R12*a1_x - R22*a1_y - R32*a1_z + R12*q_x + R22*q_y + R32*q_z;
F(20) = len1/2 + R12*a1_x + R22*a1_y + R32*a1_z - R12*q_x - R22*q_y - R32*q_z;
F(21) = wid1/2 - R13*a1_x - R23*a1_y - R33*a1_z + R13*q_x + R23*q_y + R33*q_z + (len1*wid2)/len2 - (2*wid2*(R12*a1_x + R22*a1_y + R32*a1_z - R12*q_x - R22*q_y - R32*q_z))/len2;
F(22) = wid1/2 + R13*a1_x + R23*a1_y + R33*a1_z - R13*q_x - R23*q_y - R33*q_z + (len1*wid2)/len2 + (2*wid2*(R12*a1_x + R22*a1_y + R32*a1_z - R12*q_x - R22*q_y - R32*q_z))/len2;
F(23) = wid1/2 + R13*a1_x + R23*a1_y + R33*a1_z - R13*q_x - R23*q_y - R33*q_z + (len1*wid2)/len2 - (2*wid2*(R12*a1_x + R22*a1_y + R32*a1_z - R12*q_x - R22*q_y - R32*q_z))/len2;
F(24) = wid1/2 - R13*a1_x - R23*a1_y - R33*a1_z + R13*q_x + R23*q_y + R33*q_z + (len1*wid2)/len2 + (2*wid2*(R12*a1_x + R22*a1_y + R32*a1_z - R12*q_x - R22*q_y - R32*q_z))/len2;
F(25) = wid1/2 + wid2 - R13*a1_x - R23*a1_y - R33*a1_z + R13*q_x + R23*q_y + R33*q_z;
F(26) = wid1/2 + wid2 + R13*a1_x + R23*a1_y + R33*a1_z - R13*q_x - R23*q_y - R33*q_z;
F(27) = - a2_z;

%% non-penetration 
F(28) = a1_z;

%% Jacobian matrix (Using chain rules)

if (jacflag)
    
    J1 = [                   -m,                    0,  0,                                  0,                                 0,                                  0,                         0,                         0,                         0, 0, 0,  0,              1,              0,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                                    0,          0,                                                                        0,                                        0,          0,                                                                        0,                                        0,          0,                                                                        0,                                       0,            0,            0,           0,           0,            0,            0,            0,           0,            0,                         0,                         0,                         0
                              0,                   -m,  0,                                  0,                                 0,                                  0,                         0,                         0,                         0, 0, 0,  0,              0,              1,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                                    0,          0,                                                                        0,                                        0,          0,                                                                        0,                                        0,          0,                                                                        0,                                       0,            0,            0,           0,           0,            0,            0,            0,           0,            0,                         0,                         0,                         0
                              0,                    0, -m,                                  0,                                 0,                                  0,                         0,                         0,                         0, 0, 0,  0,              0,              0,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                                    1,          0,                                                                        0,                                        0,          0,                                                                        0,                                        0,          0,                                                                        0,                                       0,            0,            0,           0,           0,            0,            0,            0,           0,            0,                         0,                         0,                         0
                              0,                    0,  0,                               -I11,                              -I12,                               -I13,                         0,                       p_n,                      -p_o, 0, 0,  0,              0,     q_z - a1_z,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                           a1_y - q_y,          0,                                                                        0,                                        0,          0,                                 -B*M_r*V_m*h*cos(phi)*sin(2*pi*fqn*time), -B*M_r*V_m*h*sin(phi)*sin(2*pi*fqn*time),          0,                                  B*M_r*V_m*h*cos(phi)*cos(2*pi*fqn*time), B*M_r*V_m*h*sin(phi)*cos(2*pi*fqn*time),   w_xo - w_x,   w_yo - w_y,  w_zo - w_z, h*w_xo*w_zo,  h*w_yo*w_zo,     h*w_zo^2, -h*w_xo*w_yo,   -h*w_yo^2, -h*w_yo*w_zo,                         0,                      -p_n,                       p_o
                              0,                    0,  0,                               -I21,                              -I22,                               -I23,                      -p_n,                         0,                       p_t, 0, 0,  0,     a1_z - q_z,              0,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                           q_x - a1_x,          0,                                  B*M_r*V_m*h*cos(phi)*sin(2*pi*fqn*time),  B*M_r*V_m*h*sin(phi)*sin(2*pi*fqn*time),          0,                                                                        0,                                        0,          0,                                                                        0,                                       0, -h*w_xo*w_zo, -h*w_yo*w_zo,   -h*w_zo^2,  w_xo - w_x,   w_yo - w_y,   w_zo - w_z,     h*w_xo^2, h*w_xo*w_yo,  h*w_xo*w_zo,                       p_n,                         0,                      -p_t
                              0,                    0,  0,                               -I31,                              -I32,                               -I33,                       p_o,                      -p_t,                         0, 0, 0,  0,     q_y - a1_y,     a1_x - q_x,              1,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                                    0,          0,                                 -B*M_r*V_m*h*cos(phi)*cos(2*pi*fqn*time), -B*M_r*V_m*h*sin(phi)*cos(2*pi*fqn*time),          0,                                                                        0,                                        0,          0,                                                                        0,                                       0,  h*w_xo*w_yo,     h*w_yo^2, h*w_yo*w_zo,   -h*w_xo^2, -h*w_xo*w_yo, -h*w_xo*w_zo,   w_xo - w_x,  w_yo - w_y,   w_zo - w_z,                      -p_o,                       p_t,                         0
           e_t^2*mu*(Van + p_n),                    0,  0,                                  0, e_t^2*mu*(Van + p_n)*(a1_z - q_z), -e_t^2*mu*(Van + p_n)*(a1_y - q_y),                         0, -e_t^2*mu*w_z*(Van + p_n),  e_t^2*mu*w_y*(Van + p_n), 0, 0,  0,            sig,              0,              0, p_t,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0, e_t^2*mu*v_x - e_t^2*mu*w_z*(a1_y - q_y) + e_t^2*mu*w_y*(a1_z - q_z),          0,                                                                        0,                                        0,          0,                                                                        0,                                        0,          0,                                                                        0,                                       0,            0,            0,           0,           0,            0,            0,            0,           0,            0,                         0,  e_t^2*mu*w_z*(Van + p_n), -e_t^2*mu*w_y*(Van + p_n)
                              0, e_o^2*mu*(Van + p_n),  0, -e_o^2*mu*(Van + p_n)*(a1_z - q_z),                                 0,  e_o^2*mu*(Van + p_n)*(a1_x - q_x),  e_o^2*mu*w_z*(Van + p_n),                         0, -e_o^2*mu*w_x*(Van + p_n), 0, 0,  0,              0,            sig,              0, p_o,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0, e_o^2*mu*v_y + e_o^2*mu*w_z*(a1_x - q_x) - e_o^2*mu*w_x*(a1_z - q_z),          0,                                                                        0,                                        0,          0,                                                                        0,                                        0,          0,                                                                        0,                                       0,            0,            0,           0,           0,            0,            0,            0,           0,            0, -e_o^2*mu*w_z*(Van + p_n),                         0,  e_o^2*mu*w_x*(Van + p_n)
                              0,                    0,  0,                                  0,                                 0,               e_r^2*mu*(Van + p_n),                         0,                         0,                         0, 0, 0,  0,              0,              0,            sig, p_r,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                         e_r^2*mu*w_z,          0,                                                                        0,                                        0,          0,                                                                        0,                                        0,          0,                                                                        0,                                       0,            0,            0,           0,           0,            0,            0,            0,           0,            0,                         0,                         0,                         0
                              0,                    0,  0,                                  0,                                 0,                                  0,                        -1,                         0,                         0, 1, 0,  0,              0,              0,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                                    0,          0,                                                                        0,                                        0,          0,                                                                        0,                                        0,          0,                                                                        0,                                       0,            0,            0,           0,           0,            0,            0,            0,           0,            0,                         0,                         0,                         0
                              0,                    0,  0,                                  0,                                 0,                                  0,                         0,                        -1,                         0, 0, 1,  0,              0,              0,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                                    0,          0,                                                                        0,                                        0,          0,                                                                        0,                                        0,          0,                                                                        0,                                       0,            0,            0,           0,           0,            0,            0,            0,           0,            0,                         0,                         0,                         0
                              0,                    0,  0,                                  0,                                 0,                                  0,                         0,                         0,                        -1, 0, 0,  1,              0,              0,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 1,                                                                    0,          0,                                                                        0,                                        0,          0,                                                                        0,                                        0,          0,                                                                        0,                                       0,            0,            0,           0,           0,            0,            0,            0,           0,            0,                         0,                         0,                         0
                              0,                    0,  0,                                  0,                                 0,                                  0,                         0,                         0,                         0, 0, 0,  0,              0,              0,              0,   0, R11, -R11, R12, -R12, R13 + (2*R12*wid2)/len2, - R13 - (2*R12*wid2)/len2, (2*R12*wid2)/len2 - R13, R13 - (2*R12*wid2)/len2, R13, -R13, 0,                                                                    0,    l1 - l2, (l3*len2 - l4*len2 + 2*l5*wid2 - 2*l6*wid2 + 2*l7*wid2 - 2*l8*wid2)/len2,             l5 - l6 - l7 + l8 + l9 - l10,          0,                                                                        0,                                        0,          0,                                                                        0,                                       0,            0,            0,           0,           0,            0,            0,            0,           0,            0,                         0,                         0,                         0
                              0,                    0,  0,                                  0,                                 0,                                  0,                         0,                         0,                         0, 0, 0,  0,              0,              0,              0,   0, R21, -R21, R22, -R22, R23 + (2*R22*wid2)/len2, - R23 - (2*R22*wid2)/len2, (2*R22*wid2)/len2 - R23, R23 - (2*R22*wid2)/len2, R23, -R23, 0,                                                                    0,          0,                                                                        0,                                        0,    l1 - l2, (l3*len2 - l4*len2 + 2*l5*wid2 - 2*l6*wid2 + 2*l7*wid2 - 2*l8*wid2)/len2,             l5 - l6 - l7 + l8 + l9 - l10,          0,                                                                        0,                                       0,            0,            0,           0,           0,            0,            0,            0,           0,            0,                         0,                         0,                         0
                              0,                    0,  0,                                  0,                                 0,                                  0,                         0,                         0,                         0, 0, 0,  0,              0,              0,              0,   0, R31, -R31, R32, -R32, R33 + (2*R32*wid2)/len2, - R33 - (2*R32*wid2)/len2, (2*R32*wid2)/len2 - R33, R33 - (2*R32*wid2)/len2, R33, -R33, 0,                                                                    0,          0,                                                                        0,                                        0,          0,                                                                        0,                                        0,    l1 - l2, (l3*len2 - l4*len2 + 2*l5*wid2 - 2*l6*wid2 + 2*l7*wid2 - 2*l8*wid2)/len2,            l5 - l6 - l7 + l8 + l9 - l10,            0,            0,           0,           0,            0,            0,            0,           0,            0,                         0,                         0,                         0
                              0,                    0,  0,                                  0,                                 0,                                  0,                         0,                         0,                         0, 0, 0,  0, -(2*p_t)/e_t^2, -(2*p_o)/e_o^2, -(2*p_r)/e_r^2,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                 mu^2*(2*Van + 2*p_n),          0,                                                                        0,                                        0,          0,                                                                        0,                                        0,          0,                                                                        0,                                       0,            0,            0,           0,           0,            0,            0,            0,           0,            0,                         0,                         0,                         0
                              0,                    0,  0,                                  0,                                 0,                                  0,                      -R11,                      -R21,                      -R31, 0, 0,  0,              0,              0,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                                    0, q_x - a1_x,                                                                        0,                                        0, q_y - a1_y,                                                                        0,                                        0, q_z - a1_z,                                                                        0,                                       0,            0,            0,           0,           0,            0,            0,            0,           0,            0,                       R11,                       R21,                       R31
                              0,                    0,  0,                                  0,                                 0,                                  0,                       R11,                       R21,                       R31, 0, 0,  0,              0,              0,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                                    0, a1_x - q_x,                                                                        0,                                        0, a1_y - q_y,                                                                        0,                                        0, a1_z - q_z,                                                                        0,                                       0,            0,            0,           0,           0,            0,            0,            0,           0,            0,                      -R11,                      -R21,                      -R31
                              0,                    0,  0,                                  0,                                 0,                                  0,                      -R12,                      -R22,                      -R32, 0, 0,  0,              0,              0,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                                    0,          0,                                                               q_x - a1_x,                                        0,          0,                                                               q_y - a1_y,                                        0,          0,                                                               q_z - a1_z,                                       0,            0,            0,           0,           0,            0,            0,            0,           0,            0,                       R12,                       R22,                       R32
                              0,                    0,  0,                                  0,                                 0,                                  0,                       R12,                       R22,                       R32, 0, 0,  0,              0,              0,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                                    0,          0,                                                               a1_x - q_x,                                        0,          0,                                                               a1_y - q_y,                                        0,          0,                                                               a1_z - q_z,                                       0,            0,            0,           0,           0,            0,            0,            0,           0,            0,                      -R12,                      -R22,                      -R32
                              0,                    0,  0,                                  0,                                 0,                                  0, - R13 - (2*R12*wid2)/len2, - R23 - (2*R22*wid2)/len2, - R33 - (2*R32*wid2)/len2, 0, 0,  0,              0,              0,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                                    0,          0,                                              -(2*wid2*(a1_x - q_x))/len2,                               q_x - a1_x,          0,                                              -(2*wid2*(a1_y - q_y))/len2,                               q_y - a1_y,          0,                                              -(2*wid2*(a1_z - q_z))/len2,                              q_z - a1_z,            0,            0,           0,           0,            0,            0,            0,           0,            0,   R13 + (2*R12*wid2)/len2,   R23 + (2*R22*wid2)/len2,   R33 + (2*R32*wid2)/len2
                              0,                    0,  0,                                  0,                                 0,                                  0,   R13 + (2*R12*wid2)/len2,   R23 + (2*R22*wid2)/len2,   R33 + (2*R32*wid2)/len2, 0, 0,  0,              0,              0,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                                    0,          0,                                               (2*wid2*(a1_x - q_x))/len2,                               a1_x - q_x,          0,                                               (2*wid2*(a1_y - q_y))/len2,                               a1_y - q_y,          0,                                               (2*wid2*(a1_z - q_z))/len2,                              a1_z - q_z,            0,            0,           0,           0,            0,            0,            0,           0,            0, - R13 - (2*R12*wid2)/len2, - R23 - (2*R22*wid2)/len2, - R33 - (2*R32*wid2)/len2
                              0,                    0,  0,                                  0,                                 0,                                  0,   R13 - (2*R12*wid2)/len2,   R23 - (2*R22*wid2)/len2,   R33 - (2*R32*wid2)/len2, 0, 0,  0,              0,              0,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                                    0,          0,                                              -(2*wid2*(a1_x - q_x))/len2,                               a1_x - q_x,          0,                                              -(2*wid2*(a1_y - q_y))/len2,                               a1_y - q_y,          0,                                              -(2*wid2*(a1_z - q_z))/len2,                              a1_z - q_z,            0,            0,           0,           0,            0,            0,            0,           0,            0,   (2*R12*wid2)/len2 - R13,   (2*R22*wid2)/len2 - R23,   (2*R32*wid2)/len2 - R33
                              0,                    0,  0,                                  0,                                 0,                                  0,   (2*R12*wid2)/len2 - R13,   (2*R22*wid2)/len2 - R23,   (2*R32*wid2)/len2 - R33, 0, 0,  0,              0,              0,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                                    0,          0,                                               (2*wid2*(a1_x - q_x))/len2,                               q_x - a1_x,          0,                                               (2*wid2*(a1_y - q_y))/len2,                               q_y - a1_y,          0,                                               (2*wid2*(a1_z - q_z))/len2,                              q_z - a1_z,            0,            0,           0,           0,            0,            0,            0,           0,            0,   R13 - (2*R12*wid2)/len2,   R23 - (2*R22*wid2)/len2,   R33 - (2*R32*wid2)/len2
                              0,                    0,  0,                                  0,                                 0,                                  0,                      -R13,                      -R23,                      -R33, 0, 0,  0,              0,              0,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                                    0,          0,                                                                        0,                               q_x - a1_x,          0,                                                                        0,                               q_y - a1_y,          0,                                                                        0,                              q_z - a1_z,            0,            0,           0,           0,            0,            0,            0,           0,            0,                       R13,                       R23,                       R33
                              0,                    0,  0,                                  0,                                 0,                                  0,                       R13,                       R23,                       R33, 0, 0,  0,              0,              0,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                                    0,          0,                                                                        0,                               a1_x - q_x,          0,                                                                        0,                               a1_y - q_y,          0,                                                                        0,                              a1_z - q_z,            0,            0,           0,           0,            0,            0,            0,           0,            0,                      -R13,                      -R23,                      -R33
                              0,                    0,  0,                                  0,                                 0,                                  0,                         0,                         0,                         0, 0, 0, -1,              0,              0,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                                    0,          0,                                                                        0,                                        0,          0,                                                                        0,                                        0,          0,                                                                        0,                                       0,            0,            0,           0,           0,            0,            0,            0,           0,            0,                         0,                         0,                         0
                              0,                    0,  0,                                  0,                                 0,                                  0,                         0,                         0,                         1, 0, 0,  0,              0,              0,              0,   0,   0,    0,   0,    0,                       0,                         0,                       0,                       0,   0,    0, 0,                                                                    0,          0,                                                                        0,                                        0,          0,                                                                        0,                                        0,          0,                                                                        0,                                       0,            0,            0,           0,           0,            0,            0,            0,           0,            0,                         0,                         0,                         0];
  
    J2 = zeros(49,32);
    J2(1:28,1:28) = eye(28); %dz/dz
    J2(29:37,end-3:end) = [2*q0,  2*q1, -2*q2, -2*q3
                          -2*q3,  2*q2,  2*q1, -2*q0
                           2*q2,  2*q3,  2*q0,  2*q1
                           2*q3,  2*q2,  2*q1,  2*q0
                           2*q0, -2*q1,  2*q2, -2*q3
                          -2*q1, -2*q0,  2*q3,  2*q2
                          -2*q2,  2*q3, -2*q0,  2*q1
                           2*q1,  2*q0,  2*q3,  2*q2
                           2*q0, -2*q1, -2*q2,  2*q3]; %dRij/d[q0,q1,q2,q3]
                       
    J2(38:46,end-3:end) = [                                                 4*I_xx*R11*q0 - 4*I_yy*R12*q3 + 4*I_zz*R13*q2,                                                 4*I_xx*R11*q1 + 4*I_yy*R12*q2 + 4*I_zz*R13*q3,                                                 4*I_yy*R12*q1 - 4*I_xx*R11*q2 + 4*I_zz*R13*q0,                                                 4*I_zz*R13*q1 - 4*I_yy*R12*q0 - 4*I_xx*R11*q3
                            2*I_xx*R11*q3 + 2*I_xx*R21*q0 + 2*I_yy*R12*q0 - 2*I_yy*R22*q3 - 2*I_zz*R13*q1 + 2*I_zz*R23*q2, 2*I_xx*R11*q2 + 2*I_xx*R21*q1 - 2*I_yy*R12*q1 + 2*I_yy*R22*q2 - 2*I_zz*R13*q0 + 2*I_zz*R23*q3, 2*I_xx*R11*q1 - 2*I_xx*R21*q2 + 2*I_yy*R12*q2 + 2*I_yy*R22*q1 + 2*I_zz*R13*q3 + 2*I_zz*R23*q0, 2*I_xx*R11*q0 - 2*I_xx*R21*q3 - 2*I_yy*R12*q3 - 2*I_yy*R22*q0 + 2*I_zz*R13*q2 + 2*I_zz*R23*q1
                            2*I_xx*R31*q0 - 2*I_xx*R11*q2 + 2*I_yy*R12*q1 - 2*I_yy*R32*q3 + 2*I_zz*R13*q0 + 2*I_zz*R33*q2, 2*I_xx*R11*q3 + 2*I_xx*R31*q1 + 2*I_yy*R12*q0 + 2*I_yy*R32*q2 - 2*I_zz*R13*q1 + 2*I_zz*R33*q3, 2*I_yy*R12*q3 - 2*I_xx*R31*q2 - 2*I_xx*R11*q0 + 2*I_yy*R32*q1 - 2*I_zz*R13*q2 + 2*I_zz*R33*q0, 2*I_xx*R11*q1 - 2*I_xx*R31*q3 + 2*I_yy*R12*q2 - 2*I_yy*R32*q0 + 2*I_zz*R13*q3 + 2*I_zz*R33*q1
                            2*I_xx*R11*q3 + 2*I_xx*R21*q0 + 2*I_yy*R12*q0 - 2*I_yy*R22*q3 - 2*I_zz*R13*q1 + 2*I_zz*R23*q2, 2*I_xx*R11*q2 + 2*I_xx*R21*q1 - 2*I_yy*R12*q1 + 2*I_yy*R22*q2 - 2*I_zz*R13*q0 + 2*I_zz*R23*q3, 2*I_xx*R11*q1 - 2*I_xx*R21*q2 + 2*I_yy*R12*q2 + 2*I_yy*R22*q1 + 2*I_zz*R13*q3 + 2*I_zz*R23*q0, 2*I_xx*R11*q0 - 2*I_xx*R21*q3 - 2*I_yy*R12*q3 - 2*I_yy*R22*q0 + 2*I_zz*R13*q2 + 2*I_zz*R23*q1
                                                                            4*I_xx*R21*q3 + 4*I_yy*R22*q0 - 4*I_zz*R23*q1,                                                 4*I_xx*R21*q2 - 4*I_yy*R22*q1 - 4*I_zz*R23*q0,                                                 4*I_xx*R21*q1 + 4*I_yy*R22*q2 + 4*I_zz*R23*q3,                                                 4*I_xx*R21*q0 - 4*I_yy*R22*q3 + 4*I_zz*R23*q2
                            2*I_xx*R31*q3 - 2*I_xx*R21*q2 + 2*I_yy*R22*q1 + 2*I_yy*R32*q0 + 2*I_zz*R23*q0 - 2*I_zz*R33*q1, 2*I_xx*R21*q3 + 2*I_xx*R31*q2 + 2*I_yy*R22*q0 - 2*I_yy*R32*q1 - 2*I_zz*R23*q1 - 2*I_zz*R33*q0, 2*I_xx*R31*q1 - 2*I_xx*R21*q0 + 2*I_yy*R22*q3 + 2*I_yy*R32*q2 - 2*I_zz*R23*q2 + 2*I_zz*R33*q3, 2*I_xx*R21*q1 + 2*I_xx*R31*q0 + 2*I_yy*R22*q2 - 2*I_yy*R32*q3 + 2*I_zz*R23*q3 + 2*I_zz*R33*q2
                            2*I_xx*R31*q0 - 2*I_xx*R11*q2 + 2*I_yy*R12*q1 - 2*I_yy*R32*q3 + 2*I_zz*R13*q0 + 2*I_zz*R33*q2, 2*I_xx*R11*q3 + 2*I_xx*R31*q1 + 2*I_yy*R12*q0 + 2*I_yy*R32*q2 - 2*I_zz*R13*q1 + 2*I_zz*R33*q3, 2*I_yy*R12*q3 - 2*I_xx*R31*q2 - 2*I_xx*R11*q0 + 2*I_yy*R32*q1 - 2*I_zz*R13*q2 + 2*I_zz*R33*q0, 2*I_xx*R11*q1 - 2*I_xx*R31*q3 + 2*I_yy*R12*q2 - 2*I_yy*R32*q0 + 2*I_zz*R13*q3 + 2*I_zz*R33*q1
                            2*I_xx*R31*q3 - 2*I_xx*R21*q2 + 2*I_yy*R22*q1 + 2*I_yy*R32*q0 + 2*I_zz*R23*q0 - 2*I_zz*R33*q1, 2*I_xx*R21*q3 + 2*I_xx*R31*q2 + 2*I_yy*R22*q0 - 2*I_yy*R32*q1 - 2*I_zz*R23*q1 - 2*I_zz*R33*q0, 2*I_xx*R31*q1 - 2*I_xx*R21*q0 + 2*I_yy*R22*q3 + 2*I_yy*R32*q2 - 2*I_zz*R23*q2 + 2*I_zz*R33*q3, 2*I_xx*R21*q1 + 2*I_xx*R31*q0 + 2*I_yy*R22*q2 - 2*I_yy*R32*q3 + 2*I_zz*R23*q3 + 2*I_zz*R33*q2
                                                                            4*I_yy*R32*q1 - 4*I_xx*R31*q2 + 4*I_zz*R33*q0,                                                 4*I_xx*R31*q3 + 4*I_yy*R32*q0 - 4*I_zz*R33*q1,                                                 4*I_yy*R32*q3 - 4*I_xx*R31*q0 - 4*I_zz*R33*q2,                                                 4*I_xx*R31*q1 + 4*I_yy*R32*q2 + 4*I_zz*R33*q3];
 
    %dIij/d[q0,q1,q2,q3]                                                                    
    J2(47:49,1:3) = h*eye(3);                   
    
    J3 = zeros(32,28);  
    J3(1:28,1:28) = eye(28);  %dz/dz
    N_J3 = (h^2*w_x^2 + h^2*w_y^2 + h^2*w_z^2 + 4)^(3/2);
    J3(29:32,4:6) = [ -(h*(q1_o*h^2*w_y^2 - q2_o*w_x*h^2*w_y + q1_o*h^2*w_z^2 - q3_o*w_x*h^2*w_z + 2*q0_o*w_x*h + 4*q1_o))/N_J3, -(h*(q2_o*h^2*w_x^2 - q1_o*w_y*h^2*w_x + q2_o*h^2*w_z^2 - q3_o*w_y*h^2*w_z + 2*q0_o*w_y*h + 4*q2_o))/N_J3, -(h*(q3_o*h^2*w_x^2 - q1_o*w_z*h^2*w_x + q3_o*h^2*w_y^2 - q2_o*w_z*h^2*w_y + 2*q0_o*w_z*h + 4*q3_o))/N_J3
                       (h*(q0_o*h^2*w_y^2 - q3_o*w_x*h^2*w_y + q0_o*h^2*w_z^2 + q2_o*w_x*h^2*w_z - 2*q1_o*w_x*h + 4*q0_o))/N_J3,  (h*(q3_o*h^2*w_x^2 - q0_o*w_y*h^2*w_x + q3_o*h^2*w_z^2 + q2_o*w_y*h^2*w_z - 2*q1_o*w_y*h + 4*q3_o))/N_J3, -(h*(q2_o*h^2*w_x^2 + q0_o*w_z*h^2*w_x + q2_o*h^2*w_y^2 + q3_o*w_z*h^2*w_y + 2*q1_o*w_z*h + 4*q2_o))/N_J3
                      -(h*(q3_o*h^2*w_y^2 + q0_o*w_x*h^2*w_y + q3_o*h^2*w_z^2 + q1_o*w_x*h^2*w_z + 2*q2_o*w_x*h + 4*q3_o))/N_J3,  (h*(q0_o*h^2*w_x^2 + q3_o*w_y*h^2*w_x + q0_o*h^2*w_z^2 - q1_o*w_y*h^2*w_z - 2*q2_o*w_y*h + 4*q0_o))/N_J3,  (h*(q1_o*h^2*w_x^2 + q3_o*w_z*h^2*w_x + q1_o*h^2*w_y^2 - q0_o*w_z*h^2*w_y - 2*q2_o*w_z*h + 4*q1_o))/N_J3
                       (h*(q2_o*h^2*w_y^2 + q1_o*w_x*h^2*w_y + q2_o*h^2*w_z^2 - q0_o*w_x*h^2*w_z - 2*q3_o*w_x*h + 4*q2_o))/N_J3, -(h*(q1_o*h^2*w_x^2 + q2_o*w_y*h^2*w_x + q1_o*h^2*w_z^2 + q0_o*w_y*h^2*w_z + 2*q3_o*w_y*h + 4*q1_o))/N_J3,  (h*(q0_o*h^2*w_x^2 - q2_o*w_z*h^2*w_x + q0_o*h^2*w_y^2 + q1_o*w_z*h^2*w_y - 2*q3_o*w_z*h + 4*q0_o))/N_J3];%d(q0 q1 q2 q3)/d(wx wy wz)
 
    J = J1*J2*J3;
    J = sparse(J);
    
end
end
