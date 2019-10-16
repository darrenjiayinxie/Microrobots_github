function A = initial_guess(A)

%% all the initial guess depends on the 1)planar sliding 2)initial oritentation with zero rotation about normal axis
unit = A.unit;
unit_mass = A.unit_mass;
infty = 1e20;
e_t = A.ellipsoid(1);
e_o = A.ellipsoid(2);
e_r = A.ellipsoid(3)*unit;



m = A.mass*unit_mass;
g = A.gravity*unit;
h = A.h;
mu = A.cof;
F_elect = A.F_elect*unit*unit_mass;

nu = A.initial_v; 
nu(1:3) = nu(1:3)*unit;
v_x = nu(1);
v_y = nu(2);
w_z = nu(3);

q_x = A.initial_q(1)*unit;
q_y = A.initial_q(2)*unit;
theta_z = A.initial_q(3);


a_y = 0; % assuming surface contact  
    
len1 = A.dim(1)*unit;
wid1 = A.dim(3)*unit;
wid2 = A.dim(4)*unit;
r1_x = A.r1_x*unit;
r = A.r*unit;

ECP1 = [q_x;a_y+wid1/2;q_x;a_y];
ECP2 = [q_x-len1/2;a_y+wid1/2-wid2/2;r1_x+r;a_y];
ECP3 = [q_x;a_y;q_x;a_y];
ECP4 = [q_x;a_y;r1_x+r;a_y];



sig1 = 0;
sig2 =0;
sig3 = v_x;
sig4 = 0;

p_n1 = 0;
p_n2 = 0;
p_n3 = (m*g*h+F_elect*h);
p_n4 = 0;


p_t1 = 0;
p_t2 = 0;
p_t3 = mu*p_n3;
p_t4 = 0;


Con_wrench = [p_t1;p_t2;p_t3;p_t4];

A.l(1:23,1) = -infty;
A.l(24:51,1) = 0;
A.u(1:51,1) = infty;
La1 =[0;0;1;0;wid1/2]; %% assuming planar contact
La2 =[0;0;0;0;0]; %% assuming planar contact
La3 =[0;0.5;0.5;0;0]; %% assuming planar contact
La4 =[0;0.5;0.5;0;0]; %% assuming planar contact


A.fun = 'mcp_funjac_microrobot_spiked_2D_2cylinders_terrain';
A.check = @mcp_funjac_microrobot_spiked_2D_2cylinders_terrain;
  
A.Z = [nu;ECP1;ECP2;ECP3;ECP4;Con_wrench;sig1;sig2;sig3;sig4;La1;La2;La3;La4;p_n1;p_n2;p_n3;p_n4]; 


end