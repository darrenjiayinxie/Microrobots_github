function [F, J, domerr] = mcp_funjac_2R_manipulator_block_sliding(z, jacflag)
%% initialize
z = z(:);
F = [];
J = [];
domerr = 0;

%% obtain value of global variables
global h;

global q_old ;
theta_1o = q_old(1);
theta_2o = q_old(2);
 



global nu_old;
w_1o =nu_old(1);
w_2o =nu_old(2);
v_xo =nu_old(3); 

 
global I_z1 I_z2 m1 m2 L1 L2 r1 r2 m H g muRB muBG eRB_t eBG_t;


global tau_1 tau_2 p_x p_y;


%% unknown variables

w_1 = z(1);
w_2 = z(2); 
v_x = z(3);

fRB_t = z(4);
fBG_t = z(5);
fBG_n = z(6);

a2RB_x = z(7); 
a2RB_y = z(8); 

lRB_k = z(9); 

fRB_n = z(10);

sigRB = z(11);  
sigBG = z(12); 



%% intermediate variables

theta1 = theta_1o + h*w_1;
theta2 = theta_2o + h*w_2;



%% simplification
alpha = I_z1 + I_z2 + m1*r1^2 +m2*(L1^2 +r2^2);
beta = m2*L1*r2;
delta = I_z2+m2*r2^2;



%% Dynamic equations
F(1) = h*(g*m2*r2*cos(theta1 + theta2) + L1*g*m2*cos(theta1) + g*m1*r1*cos(theta1)) - h*tau_1 - h*(beta*w_2*sin(theta2)*(w_1 + w_2) + beta*w_1*w_2*sin(theta2)) + (alpha + 2*beta*cos(theta2))*(w_1 - w_1o) + (delta + beta*cos(theta2))*(w_2 - w_2o) - fRB_n*h*(L2*cos(theta1 + theta2) + L1*cos(theta1)) + fRB_t*h*(L2*sin(theta1 + theta2) + L1*sin(theta1));
F(2) = (delta + beta*cos(theta2))*(w_1 - w_1o) - h*tau_2 + delta*(w_2 - w_2o) - L2*fRB_n*h*cos(theta1 + theta2) + L2*fRB_t*h*sin(theta1 + theta2) + beta*h*w_1^2*sin(theta2) + g*h*m2*r2*cos(theta1 + theta2);

F(3) = fRB_t*h - fBG_t*h - p_x + m*(v_x - v_xo);
F(4) = fRB_n*h - fBG_n*h - p_y + g*h*m;

%% Friction model without complementarity equation
F(5) = - fRB_n*h*muRB*(v_x + w_1*(L2*sin(theta1 + theta2) + L1*sin(theta1)) + L2*w_2*sin(theta1 + theta2))*eRB_t^2 + fRB_t*h*sigRB;
F(6) = fBG_n*h*muBG*v_x*eBG_t^2 + fBG_t*h*sigBG;

%% contact constraints
F(7) = L2*cos(theta1 + theta2) - a2RB_x + L1*cos(theta1);
F(8) = L2*sin(theta1 + theta2) - lRB_k - a2RB_y + L1*sin(theta1);

F(9) =  H - a2RB_y;
F(10) = L2*sin(theta1 + theta2) - H + L1*sin(theta1);

%% Friction model's complementarity equation
F(11) = fRB_n^2*muRB^2 - fRB_t^2/eRB_t^2;
F(12) = fBG_n^2*muBG^2 - fBG_t^2/eBG_t^2;

 

if (jacflag)
    %% J1
     J1 = [            alpha + 2*beta*cos(theta2) - 2*beta*h*w_2*sin(theta2), delta - h*(beta*w_1*sin(theta2) + beta*w_2*sin(theta2) + beta*sin(theta2)*(w_1 + w_2)) + beta*cos(theta2),                     0, h*(L2*sin(theta1 + theta2) + L1*sin(theta1)),                  0,                  0,  0,  0,  0,                                                        -h*(L2*cos(theta1 + theta2) + L1*cos(theta1)),       0,       0, fRB_t*h*(L2*cos(theta1 + theta2) + L1*cos(theta1)) - h*(g*m2*r2*sin(theta1 + theta2) + L1*g*m2*sin(theta1) + g*m1*r1*sin(theta1)) + fRB_n*h*(L2*sin(theta1 + theta2) + L1*sin(theta1)), L2*fRB_t*h*cos(theta1 + theta2) - 2*beta*sin(theta2)*(w_1 - w_1o) - beta*sin(theta2)*(w_2 - w_2o) - h*(beta*w_2*cos(theta2)*(w_1 + w_2) + beta*w_1*w_2*cos(theta2)) + L2*fRB_n*h*sin(theta1 + theta2) - g*h*m2*r2*sin(theta1 + theta2)
                         delta + beta*cos(theta2) + 2*beta*h*w_1*sin(theta2),                                                                                                     delta,                     0,                    L2*h*sin(theta1 + theta2),                  0,                  0,  0,  0,  0,                                                                           -L2*h*cos(theta1 + theta2),       0,       0,                                                                                     L2*fRB_t*h*cos(theta1 + theta2) + L2*fRB_n*h*sin(theta1 + theta2) - g*h*m2*r2*sin(theta1 + theta2),                                                                          L2*fRB_t*h*cos(theta1 + theta2) - beta*sin(theta2)*(w_1 - w_1o) + L2*fRB_n*h*sin(theta1 + theta2) + beta*h*w_1^2*cos(theta2) - g*h*m2*r2*sin(theta1 + theta2)
                                                                           0,                                                                                                         0,                     m,                                            h,                 -h,                  0,  0,  0,  0,                                                                                                    0,       0,       0,                                                                                                                                                                                      0,                                                                                                                                                                                                                                      0
	                                                                       0,                                                                                                         0,                     0,                                            0,                  0,                 -h,  0,  0,  0,                                                                                                    h,       0,       0,                                                                                                                                                                                      0,                                                                                                                                                                                                                                      0
            -eRB_t^2*fRB_n*h*muRB*(L2*sin(theta1 + theta2) + L1*sin(theta1)),                                                             -L2*eRB_t^2*fRB_n*h*muRB*sin(theta1 + theta2), -eRB_t^2*fRB_n*h*muRB,                                      h*sigRB,                  0,                  0,  0,  0,  0, -eRB_t^2*h*muRB*(v_x + w_1*(L2*sin(theta1 + theta2) + L1*sin(theta1)) + L2*w_2*sin(theta1 + theta2)), fRB_t*h,       0,                                                                                   -eRB_t^2*fRB_n*h*muRB*(w_1*(L2*cos(theta1 + theta2) + L1*cos(theta1)) + L2*w_2*cos(theta1 + theta2)),                                                                                                                                                      -eRB_t^2*fRB_n*h*muRB*(L2*w_1*cos(theta1 + theta2) + L2*w_2*cos(theta1 + theta2))
                                                                           0,                                                                                                         0,  eBG_t^2*fBG_n*h*muBG,                                            0,            h*sigBG, eBG_t^2*h*muBG*v_x,  0,  0,  0,                                                                                                    0,       0, fBG_t*h,                                                                                                                                                                                      0,                                                                                                                                                                                                                                      0
                                                                           0,                                                                                                         0,                     0,                                            0,                  0,                  0, -1,  0,  0,                                                                                                    0,       0,       0,                                                                                                                                             - L2*sin(theta1 + theta2) - L1*sin(theta1),                                                                                                                                                                                                               -L2*sin(theta1 + theta2)
                                                                           0,                                                                                                         0,                     0,                                            0,                  0,                  0,  0, -1, -1,                                                                                                    0,       0,       0,                                                                                                                                               L2*cos(theta1 + theta2) + L1*cos(theta1),                                                                                                                                                                                                                L2*cos(theta1 + theta2)
                                                                           0,                                                                                                         0,                     0,                                            0,                  0,                  0,  0, -1,  0,                                                                                                    0,       0,       0,                                                                                                                                                                                      0,                                                                                                                                                                                                                                      0
                                                                           0,                                                                                                         0,                     0,                                            0,                  0,                  0,  0,  0,  0,                                                                                                    0,       0,       0,                                                                                                                                               L2*cos(theta1 + theta2) + L1*cos(theta1),                                                                                                                                                                                                                L2*cos(theta1 + theta2)
                                                                           0,                                                                                                         0,                     0,                           -(2*fRB_t)/eRB_t^2,                  0,                  0,  0,  0,  0,                                                                                       2*fRB_n*muRB^2,       0,       0,                                                                                                                                                                                      0,                                                                                                                                                                                                                                      0
                                                                           0,                                                                                                         0,                     0,                                            0, -(2*fBG_t)/eBG_t^2,     2*fBG_n*muBG^2,  0,  0,  0,                                                                                                    0,       0,       0,                                                                                                                                                                                      0,                                                                                                                                                                                                                                      0];
                                                                       
     J2 = zeros(14,12);
     J2(1:12,1:12) = eye(12);
     J2(13:14,1:2) = h*eye(2);
     J = J1*J2;
     J = sparse(J);
end








end