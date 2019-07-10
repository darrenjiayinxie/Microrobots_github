Van = 3.430826822607249e-06*(A2.unit*A2.unit_mass*A2.h);
sig = A2.z(16,38);
p_o = A2.z(14,38);
e_t = 1;
e_o = 1;
e_r = 0.1;
mu = A2.cof;
v_y = A2.z(2,38);
theta = A2.theta;
p_n = A2.z(28,38);
v_z = A2.z(3,38);
w_x = A2.z(4,38);
w_y = A2.z(5,38);
w_z = A2.z(6,38);
a1_z = A2.z(9,38);
q_z = A2.q(3,38);
a1_y = A2.z(8,38);
q_y = A2.q(2,38);
a1_x = A2.z(7,38);
q_x = A2.q(1,38);


F(8) = p_o*sig + e_o^2*mu*v_y*cos(theta)*(Van + p_n) - e_o^2*mu*v_z*sin(theta)*(Van + p_n) - e_o^2*mu*w_x*(cos(theta)*(a1_z - q_z) + sin(theta)*(a1_y - q_y))*(Van + p_n) + e_o^2*mu*w_z*cos(theta)*(Van + p_n)*(a1_x - q_x) + e_o^2*mu*w_y*sin(theta)*(Van + p_n)*(a1_x - q_x);