function [F, J, domerr] = mcp_funjac_planar_pushing(z, jacflag)
%% initialize
z = z(:);
F = [];
J = [];
domerr = 0;
%% obtain value of global variables
global h;

global q_old;
q_xo = q_old(1);
q_yo = q_old(2);
theta_zo = q_old(3);
q_xso = q_old(4);
q_yso = q_old(5);
theta_zso = q_old(6);

global nu_old;
v_xo =nu_old(1);
v_yo =nu_old(2);
w_zo =nu_old(3);
v_xso =nu_old(4);
v_yso =nu_old(5);
w_zso =nu_old(6);

global I_z m I_zs m_s;

global mu mu_s e_t e_o e_r;

global q_z q_zs p_ns;
 
global p_x p_y p_zt; % applied impulse for the pusher

global len_s; % dimension of the cuboid

global len wid theta_x; %dimension of the cylinder

%% unknown variables
v_x = z(1);
v_y = z(2);
w_z = z(3);
v_xs = z(4);
v_ys = z(5);
w_zs = z(6);
p_t = z(7);
p_o = z(8);
p_t_s = z(9);
p_o_s = z(10);
p_r_s = z(11);
a1_x = z(12);
a1_y = z(13);
a1_z = z(14);
a2_x = z(15);
a2_y = z(16);
a2_z = z(17);
sig = z(18);
sig_s = z(19);
l1 = z(20);
l2 = z(21);
l3 = z(22);
l4 = z(23);
l5 = z(24);
l6 = z(25);
l7 = z(26);
p_n = z(27);

%% intermediate variables
q_x = q_xo+h*v_x;
q_y = q_yo+h*v_y;
theta_z = theta_zo + h*w_z;

q_xs = q_xso+h*v_xs;
q_ys = q_yso+h*v_ys;
theta_zs = theta_zso+h*w_zs;

%% variables for simplifcation
a2_bx = a2_x*cos(theta_z) + a2_y*sin(theta_z) - q_x*cos(theta_z) - q_y*sin(theta_z);
a2_by = a2_y*cos(theta_z) - a2_x*sin(theta_z) - q_y*cos(theta_z) + q_x*sin(theta_z);
a1_bsx = cos(theta_zs)*(a1_x - q_xs) + sin(theta_zs)*(a1_y - q_ys);
a1_bsy = a1_y*cos(theta_zs) - a1_x*sin(theta_zs) - q_ys*cos(theta_zs) + q_xs*sin(theta_zs);
a2_bsx = cos(theta_zs)*(a2_x - q_xs) + sin(theta_zs)*(a2_y - q_ys);
a2_bsy = a2_y*cos(theta_zs) - a2_x*sin(theta_zs) - q_ys*cos(theta_zs) + q_xs*sin(theta_zs);
p_xs = -(p_t*cos(theta_zs) + p_n*sin(theta_zs));
p_ys = p_n*cos(theta_zs) - p_t*sin(theta_zs);
r_n = cos(theta_zs)*(a2_x - q_x) + sin(theta_zs)*(a2_y - q_y);
r_t = cos(theta_zs)*(a2_y - q_y) - sin(theta_zs)*(a2_x - q_x);
v_t = cos(theta_zs)*(v_x - v_xs - w_z*(a2_y - q_y) + w_zs*(a1_y - q_ys)) + sin(theta_zs)*(v_y - v_ys + w_z*(a2_x - q_x) - w_zs*(a1_x - q_xs));
v_o = cos(theta_zs)*(v_y - v_ys + w_z*(a2_x - q_x) - w_zs*(a1_x - q_xs)) - sin(theta_zs)*(v_x - v_xs - w_z*(a2_y - q_y) + w_zs*(a1_y - q_ys));
p_zts = (p_ys*(a1_x-q_xs) - p_xs*(a1_y-q_ys) + ((p_xs*(p_o_s*q_zs - p_ys*(a1_z - q_zs))*(a1_z - q_zs))/p_ns - (p_ys*(p_t_s*q_zs - p_xs*(a1_z - q_zs))*(a1_z - q_zs))/p_ns)/q_zs);

A1 = - l6*cos(theta_z) + l7*cos(theta_z) - l2*cos(theta_x)*sin(theta_z) + l4*cos(theta_x)*sin(theta_z) + l3*sin(theta_x)*sin(theta_z) - l5*sin(theta_x)*sin(theta_z);
A2 =   l7*sin(theta_z) - l6*sin(theta_z) + l2*cos(theta_x)*cos(theta_z) - l4*cos(theta_x)*cos(theta_z) - l3*cos(theta_z)*sin(theta_x) + l5*cos(theta_z)*sin(theta_x);

B1 = ((a1_y*p_ns*cos(theta_zs) - a1_z*p_o_s*cos(theta_zs) - a1_x*p_ns*sin(theta_zs) + a1_z*p_t_s*sin(theta_zs) - p_ns*q_ys*cos(theta_zs) + p_o_s*q_zs*cos(theta_zs) + p_ns*q_xs*sin(theta_zs) - p_t_s*q_zs*sin(theta_zs))/p_ns);
B2 = (-(p_n*p_t_s*cos(theta_zs) + p_o_s*p_t*cos(theta_zs) + p_n*p_o_s*sin(theta_zs) - p_t*p_t_s*sin(theta_zs))/p_ns);
B3 = (a1_x*p_ns*cos(theta_zs) - a1_z*p_t_s*cos(theta_zs) + a1_y*p_ns*sin(theta_zs) - a1_z*p_o_s*sin(theta_zs) - p_ns*q_xs*cos(theta_zs) + p_t_s*q_zs*cos(theta_zs) - p_ns*q_ys*sin(theta_zs) + p_o_s*q_zs*sin(theta_zs))/p_ns;
B4 = p_xs*(a1_x - q_xs) + p_ys*(a1_y - q_ys) - ((p_xs^2*(a1_z - q_zs)^2)/p_ns + (p_ys^2*(a1_z - q_zs)^2)/p_ns + (p_ys*(p_o_s*q_zs - p_ys*(a1_z - q_zs))*(a1_z - q_zs))/p_ns + (p_xs*(p_t_s*q_zs - p_xs*(a1_z - q_zs))*(a1_z - q_zs))/p_ns)/q_zs;

%% functions 
F(1) = p_x - p_xs - m*(v_x - v_xo);
F(2) = p_y - p_ys - m*(v_y - v_yo);
F(3) = p_zt - I_z*(w_z - w_zo) - p_ys*(a2_x - q_x) + p_xs*(a2_y - q_y);

F(4) = p_t_s + p_xs - m_s*(v_xs - v_xso);
F(5) = p_o_s + p_ys - m_s*(v_ys - v_yso);
F(6) = p_r_s - I_zs*(w_zs - w_zso) + p_zts;

F(7) = mu*p_n*v_t*e_t^2 + p_t*sig;
F(8) = p_o*sig;

F(9) = mu_s*p_ns*(v_xs + (w_zs*(p_o_s*q_zs - p_ys*(a1_z - q_zs)))/p_ns)*e_t^2 + p_t_s*sig_s;
F(10) = mu_s*p_ns*(v_ys - (w_zs*(p_t_s*q_zs - p_xs*(a1_z - q_zs)))/p_ns)*e_o^2 + p_o_s*sig_s;
F(11) = mu_s*p_ns*w_zs*e_r^2 + p_r_s*sig_s;

F(12) = a1_x - a2_x + l1*sin(theta_zs);
F(13) = a1_y - a2_y - l1*cos(theta_zs);
F(14) = a1_z - a2_z;

F(15) = A1 + sin(theta_zs);
F(16) = A2 - cos(theta_zs);
F(17) = l3*cos(theta_x) - l5*cos(theta_x) + l2*sin(theta_x) - l4*sin(theta_x);

F(18) = mu^2*p_n^2 - p_t^2/e_t^2 - p_o^2/e_o^2;
F(19) = mu_s^2*p_ns^2 - p_r_s^2/e_r^2 - p_t_s^2/e_t^2 - p_o_s^2/e_o^2;

F(20) = a1_bsy + len_s/2;

F(21) = len/2 - a2_by*cos(theta_x) - sin(theta_x)*(a2_z - q_z);
F(22) = len/2 + a2_by*sin(theta_x) - cos(theta_x)*(a2_z - q_z);
F(23) = len/2 + a2_by*cos(theta_x) + sin(theta_x)*(a2_z - q_z);
F(24) = len/2 - a2_by*sin(theta_x) + cos(theta_x)*(a2_z - q_z);
F(25) = a2_bx + wid/2;
F(26) = wid/2 - a2_bx;

F(27) = - a2_bsy - len_s/2;



if(jacflag)
    J1 = [                         -m,                          0,                 0,                           0,                           0,                                             0,                                cos(theta_zs),              0,                          0,                         0,                0,                                0,                               0,                     0,                              0,                               0,             0,   0,     0,              0,                          0,                          0,                          0,                          0,             0,            0,                                sin(theta_zs),                               0,                              0,                   0,                               0,                                0,                                p_ys
                                    0,                         -m,                 0,                           0,                           0,                                             0,                                sin(theta_zs),              0,                          0,                         0,                0,                                0,                               0,                     0,                              0,                               0,             0,   0,     0,              0,                          0,                          0,                          0,                          0,             0,            0,                               -cos(theta_zs),                               0,                              0,                   0,                               0,                                0,                               -p_xs
                                    0,                          0,              -I_z,                           0,                           0,                                             0,                                         -r_t,              0,                          0,                         0,                0,                                0,                               0,                     0,                          -p_ys,                            p_xs,             0,   0,     0,              0,                          0,                          0,                          0,                          0,             0,            0,                                         -r_n,                            p_ys,                          -p_xs,                   0,                               0,                                0,                   p_t*r_n - p_n*r_t
                                    0,                          0,                 0,                        -m_s,                           0,                                             0,                               -cos(theta_zs),              0,                          1,                         0,                0,                                0,                               0,                     0,                              0,                               0,             0,   0,     0,              0,                          0,                          0,                          0,                          0,             0,            0,                               -sin(theta_zs),                               0,                              0,                   0,                               0,                                0,                               -p_ys
                                    0,                          0,                 0,                           0,                        -m_s,                                             0,                               -sin(theta_zs),              0,                          0,                         1,                0,                                0,                               0,                     0,                              0,                               0,             0,   0,     0,              0,                          0,                          0,                          0,                          0,             0,            0,                                cos(theta_zs),                               0,                              0,                   0,                               0,                                0,                                p_xs
                                    0,                          0,                 0,                           0,                           0,                                         -I_zs,                                           B1,              0, -(p_ys*(a1_z - q_zs))/p_ns, (p_xs*(a1_z - q_zs))/p_ns,                1,                             p_ys,                           -p_xs,                    B2,                              0,                               0,             0,   0,     0,              0,                          0,                          0,                          0,                          0,             0,            0,                                           B3,                               0,                              0,                   0,                           -p_ys,                             p_xs,                                  B4
           e_t^2*mu*p_n*cos(theta_zs), e_t^2*mu*p_n*sin(theta_zs), -e_t^2*mu*p_n*r_t, -e_t^2*mu*p_n*cos(theta_zs), -e_t^2*mu*p_n*sin(theta_zs),                           a1_bsy*e_t^2*mu*p_n,                                          sig,              0,                          0,                         0,                0, -e_t^2*mu*p_n*w_zs*sin(theta_zs), e_t^2*mu*p_n*w_zs*cos(theta_zs),                     0, e_t^2*mu*p_n*w_z*sin(theta_zs), -e_t^2*mu*p_n*w_z*cos(theta_zs),             0, p_t,     0,              0,                          0,                          0,                          0,                          0,             0,            0,                                 e_t^2*mu*v_t, -e_t^2*mu*p_n*w_z*sin(theta_zs), e_t^2*mu*p_n*w_z*cos(theta_zs),                   0, e_t^2*mu*p_n*w_zs*sin(theta_zs), -e_t^2*mu*p_n*w_zs*cos(theta_zs),                   -e_t^2*mu*p_n*v_o
                                    0,                          0,                 0,                           0,                           0,                                             0,                                            0,            sig,                          0,                         0,                0,                                0,                               0,                     0,                              0,                               0,             0, p_o,     0,              0,                          0,                          0,                          0,                          0,             0,            0,                                            0,                               0,                              0,                   0,                               0,                                0,                                   0
                                    0,                          0,                 0,             e_t^2*mu_s*p_ns,                           0,  e_t^2*mu_s*(p_o_s*q_zs - p_ys*(a1_z - q_zs)),  e_t^2*mu_s*w_zs*sin(theta_zs)*(a1_z - q_zs),              0,                      sig_s,      e_t^2*mu_s*q_zs*w_zs,                0,                                0,                               0, -e_t^2*mu_s*p_ys*w_zs,                              0,                               0,             0,   0, p_t_s,              0,                          0,                          0,                          0,                          0,             0,            0, -e_t^2*mu_s*w_zs*cos(theta_zs)*(a1_z - q_zs),                               0,                              0,                   0,                               0,                                0, -e_t^2*mu_s*p_xs*w_zs*(a1_z - q_zs)
                                    0,                          0,                 0,                           0,             e_o^2*mu_s*p_ns, -e_o^2*mu_s*(p_t_s*q_zs - p_xs*(a1_z - q_zs)), -e_o^2*mu_s*w_zs*cos(theta_zs)*(a1_z - q_zs),              0,      -e_o^2*mu_s*q_zs*w_zs,                     sig_s,                0,                                0,                               0,  e_o^2*mu_s*p_xs*w_zs,                              0,                               0,             0,   0, p_o_s,              0,                          0,                          0,                          0,                          0,             0,            0, -e_o^2*mu_s*w_zs*sin(theta_zs)*(a1_z - q_zs),                               0,                              0,                   0,                               0,                                0, -e_o^2*mu_s*p_ys*w_zs*(a1_z - q_zs)
                                    0,                          0,                 0,                           0,                           0,                               e_r^2*mu_s*p_ns,                                            0,              0,                          0,                         0,            sig_s,                                0,                               0,                     0,                              0,                               0,             0,   0, p_r_s,              0,                          0,                          0,                          0,                          0,             0,            0,                                            0,                               0,                              0,                   0,                               0,                                0,                                   0
                                    0,                          0,                 0,                           0,                           0,                                             0,                                            0,              0,                          0,                         0,                0,                                1,                               0,                     0,                             -1,                               0,             0,   0,     0,  sin(theta_zs),                          0,                          0,                          0,                          0,             0,            0,                                            0,                               0,                              0,                   0,                               0,                                0,                    l1*cos(theta_zs)
                                    0,                          0,                 0,                           0,                           0,                                             0,                                            0,              0,                          0,                         0,                0,                                0,                               1,                     0,                              0,                              -1,             0,   0,     0, -cos(theta_zs),                          0,                          0,                          0,                          0,             0,            0,                                            0,                               0,                              0,                   0,                               0,                                0,                    l1*sin(theta_zs)
                                    0,                          0,                 0,                           0,                           0,                                             0,                                            0,              0,                          0,                         0,                0,                                0,                               0,                     1,                              0,                               0,            -1,   0,     0,              0,                          0,                          0,                          0,                          0,             0,            0,                                            0,                               0,                              0,                   0,                               0,                                0,                                   0
                                    0,                          0,                 0,                           0,                           0,                                             0,                                            0,              0,                          0,                         0,                0,                                0,                               0,                     0,                              0,                               0,             0,   0,     0,              0, -cos(theta_x)*sin(theta_z),  sin(theta_x)*sin(theta_z),  cos(theta_x)*sin(theta_z), -sin(theta_x)*sin(theta_z), -cos(theta_z), cos(theta_z),                                            0,                               0,                              0,                 -A2,                               0,                                0,                       cos(theta_zs)
                                    0,                          0,                 0,                           0,                           0,                                             0,                                            0,              0,                          0,                         0,                0,                                0,                               0,                     0,                              0,                               0,             0,   0,     0,              0,  cos(theta_x)*cos(theta_z), -cos(theta_z)*sin(theta_x), -cos(theta_x)*cos(theta_z),  cos(theta_z)*sin(theta_x), -sin(theta_z), sin(theta_z),                                            0,                               0,                              0,                  A1,                               0,                                0,                       sin(theta_zs)
                                    0,                          0,                 0,                           0,                           0,                                             0,                                            0,              0,                          0,                         0,                0,                                0,                               0,                     0,                              0,                               0,             0,   0,     0,              0,               sin(theta_x),               cos(theta_x),              -sin(theta_x),              -cos(theta_x),             0,            0,                                            0,                               0,                              0,                   0,                               0,                                0,                                   0
                                    0,                          0,                 0,                           0,                           0,                                             0,                               -(2*p_t)/e_t^2, -(2*p_o)/e_o^2,                          0,                         0,                0,                                0,                               0,                     0,                              0,                               0,             0,   0,     0,              0,                          0,                          0,                          0,                          0,             0,            0,                                   2*mu^2*p_n,                               0,                              0,                   0,                               0,                                0,                                   0
                                    0,                          0,                 0,                           0,                           0,                                             0,                                            0,              0,           -(2*p_t_s)/e_t^2,          -(2*p_o_s)/e_o^2, -(2*p_r_s)/e_r^2,                                0,                               0,                     0,                              0,                               0,             0,   0,     0,              0,                          0,                          0,                          0,                          0,             0,            0,                                            0,                               0,                              0,                   0,                               0,                                0,                                   0
                                    0,                          0,                 0,                           0,                           0,                                             0,                                            0,              0,                          0,                         0,                0,                   -sin(theta_zs),                   cos(theta_zs),                     0,                              0,                               0,             0,   0,     0,              0,                          0,                          0,                          0,                          0,             0,            0,                                            0,                               0,                              0,                   0,                   sin(theta_zs),                   -cos(theta_zs),                             -a2_bsx
                                    0,                          0,                 0,                           0,                           0,                                             0,                                            0,              0,                          0,                         0,                0,                                0,                               0,                     0,      cos(theta_x)*sin(theta_z),      -cos(theta_x)*cos(theta_z), -sin(theta_x),   0,     0,              0,                          0,                          0,                          0,                          0,             0,            0,                                            0,      -cos(theta_x)*sin(theta_z),      cos(theta_x)*cos(theta_z),  a2_bx*cos(theta_x),                               0,                                0,                                   0
                                    0,                          0,                 0,                           0,                           0,                                             0,                                            0,              0,                          0,                         0,                0,                                0,                               0,                     0,     -sin(theta_x)*sin(theta_z),       cos(theta_z)*sin(theta_x), -cos(theta_x),   0,     0,              0,                          0,                          0,                          0,                          0,             0,            0,                                            0,       sin(theta_x)*sin(theta_z),     -cos(theta_z)*sin(theta_x), -a2_bx*sin(theta_x),                               0,                                0,                                   0
                                    0,                          0,                 0,                           0,                           0,                                             0,                                            0,              0,                          0,                         0,                0,                                0,                               0,                     0,     -cos(theta_x)*sin(theta_z),       cos(theta_x)*cos(theta_z),  sin(theta_x),   0,     0,              0,                          0,                          0,                          0,                          0,             0,            0,                                            0,       cos(theta_x)*sin(theta_z),     -cos(theta_x)*cos(theta_z), -a2_bx*cos(theta_x),                               0,                                0,                                   0
                                    0,                          0,                 0,                           0,                           0,                                             0,                                            0,              0,                          0,                         0,                0,                                0,                               0,                     0,      sin(theta_x)*sin(theta_z),      -cos(theta_z)*sin(theta_x),  cos(theta_x),   0,     0,              0,                          0,                          0,                          0,                          0,             0,            0,                                            0,      -sin(theta_x)*sin(theta_z),      cos(theta_z)*sin(theta_x),  a2_bx*sin(theta_x),                               0,                                0,                                   0
                                    0,                          0,                 0,                           0,                           0,                                             0,                                            0,              0,                          0,                         0,                0,                                0,                               0,                     0,                   cos(theta_z),                    sin(theta_z),             0,   0,     0,              0,                          0,                          0,                          0,                          0,             0,            0,                                            0,                   -cos(theta_z),                  -sin(theta_z),               a2_by,                               0,                                0,                                   0
                                    0,                          0,                 0,                           0,                           0,                                             0,                                            0,              0,                          0,                         0,                0,                                0,                               0,                     0,                  -cos(theta_z),                   -sin(theta_z),             0,   0,     0,              0,                          0,                          0,                          0,                          0,             0,            0,                                            0,                    cos(theta_z),                   sin(theta_z),              -a2_by,                               0,                                0,                                   0
                                    0,                          0,                 0,                           0,                           0,                                             0,                                            0,              0,                          0,                         0,                0,                                0,                               0,                     0,                  sin(theta_zs),                  -cos(theta_zs),             0,   0,     0,              0,                          0,                          0,                          0,                          0,             0,            0,                                            0,                               0,                              0,                   0,                  -sin(theta_zs),                    cos(theta_zs),                              a2_bsx];
    J2 = zeros(32,27);        
    J2(1:27,1:27) = eye(27);
    J2(28:33,1:6) = h*eye(6); 
    J = J1*J2;
    J = sparse(J);
end
    
    
    
    
    
end