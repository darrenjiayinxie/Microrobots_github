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
w_x = nu(4);
w_y = nu(5);
w_z = nu(6);

q_x = A.initial_q(1)*unit;
q_y = A.initial_q(2)*unit;
q_z = A.initial_q(3)*unit;


a_z = 0; % assuming surface contact  
    
    

   

len = A.dim(1)*unit;

ECP1 = [q_x;q_y-len/2;a_z;q_x;q_y-len/2;a_z];

ECP2 = [q_x;q_y-len/2;a_z;q_x;(A.r1_y+A.rc_1)*unit;a_z];
ECP3 = [q_x;q_y-len/2;a_z;q_x;(A.r2_y+A.rc_2)*unit;a_z];

v_t1 = v_x - w_z*(ECP1(2) - q_y) + w_y*(ECP1(3) - q_z);
v_o1 = v_y + w_z*(ECP1(1) - q_x) - w_x*(ECP1(3) - q_z);
v_r1 = w_z;

sig1 = sqrt(e_t^2*v_t1^2 + e_o^2*v_o1^2 + e_r^2*v_r1^2);
sig2 =0;
sig3 = 0;

p_n1 = (m*g*h+F_elect*h);
p_n2 = 0;
p_n3 = 0;
if sig1 == 0
    p_t1 = 0;
    p_o1 = 0;
    p_r1 = 0;
else
    p_t1 = -e_t^2*mu*p_n1*v_t1/sig1;
    p_o1 = -e_o^2*mu*p_n1*v_o1/sig1;
    p_r1 = -e_r^2*mu*p_n1*v_r1/sig1;
end




Con_wrench1 = [p_t1;p_o1;p_r1];

Con_wrench2 = [0;0;0];
Con_wrench3 = [0;0;0];

A.l(1:33,1) = -infty;
A.l(34:60,1) = 0;
A.u(1:60,1) = infty;
La1 =[0;0;0;0;0;1;0]; %% assuming planar contact
La2 =[0;0;0;0.4;0;0;0.3750]; %% assuming planar contact
La3 =[0;0;0;0.4;0;0;1.6250]; %% assuming planar contact

A.fun = 'mcp_funjac_microrobot_cuboid_2cylinders_terrain';
A.check = @mcp_funjac_microrobot_cuboid_2cylinders_terrain;
  
A.Z = [nu;ECP1;ECP2;ECP3;Con_wrench1;Con_wrench2;Con_wrench3;sig1;sig2;sig3;La1;p_n1;La2;p_n2;La3;p_n3]; 


end