function [F, J, domerr] = mcp_funjac_planar_sliding_convex_patches_2(z, jacflag)

%% initialize
z = z(:);
F = [];
J = [];
domerr = 0;

global v_xo v_yo w_zo;

global I_z m q_z p_n e_t e_o e_r mu;

global p_x p_y p_zt;
p_xt = 0;
p_yt = 0;

F(1) = p1_t + p2_t  + p_x - m*(v_x - v_xo);
F(2) = p1_o + p2_o + p_y - m*(v_y - v_yo);
F(3) = p1_n + p2_n - g*h*m;

F(4) = p_xt - p1_o*(0 - q_z) - p2_o*(0 - q_z) + p1_n*(a11_y - q_y) + p2_n*(a21_y - q_y);
F(5) = p_yt - p1_n*(a11_x - q_x) - p2_n*(a21_x - q_x) + p1_t*(0 - q_z) + p2_t*(0 - q_z);
F(6) = p1_r + p2_r + p_zt + p1_o*(a11_x - q_x) + p2_o*(a21_x - q_x) - p1_t*(a11_y - q_y) - p2_t*(a21_y - q_y) - I33*(w_z - w_zo);
