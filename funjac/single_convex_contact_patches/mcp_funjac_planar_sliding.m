function [F, J, domerr] = mcp_funjac_planar_sliding(z, jacflag)
%% initialize
z = z(:);
F = [];
J = [];
domerr = 0;

global v_xo v_yo w_zo;

global I_z m q_z p_n e_t e_o e_r mu;

global p_x p_y p_zt;

p_t = z(1);
p_o = z(2);
p_r = z(3);
sig = z(4);

p_xt = 0;
p_yt = 0;

v_x = (p_t+p_x)/m +v_xo;
v_y = (p_o+p_y)/m +v_yo;
w_z = (p_r+p_zt)/I_z +w_zo;

r_x = -p_t*q_z/p_n;
r_y = -p_o*q_z/p_n;

v_t = v_x -w_z*r_y;
v_o = v_y +w_z*r_x;
v_r = w_z;

F(1) = e_t^2*mu*p_n*v_t +p_t*sig;
F(2) = e_o^2*mu*p_n*v_o +p_o*sig;
F(3) = e_r^2*mu*p_n*v_r +p_r*sig;
F(4) = mu^2*p_n^2-p_t^2/e_t^2-p_o^2/e_o^2-p_r^2/e_r^2;

if(jacflag)
    J = [                (mu*p_n*e_t^2)/m + sig, e_t^2*mu*q_z*(w_zo + (p_r + p_zt)/I_z), (e_t^2*mu*(p_xt + p_o*q_z))/I_z, p_t
          -e_o^2*mu*q_z*(w_zo + (p_r + p_zt)/I_z),               (mu*p_n*e_o^2)/m + sig, (e_o^2*mu*(p_yt - p_t*q_z))/I_z, p_o
                                                0,                                      0,      (mu*p_n*e_r^2)/I_z + sig, p_r
                                   -(2*p_t)/e_t^2,                         -(2*p_o)/e_o^2,                  -(2*p_r)/e_r^2,   0];

    J = sparse(J);
end
    
    
    
    
    
end