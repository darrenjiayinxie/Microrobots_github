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

if A.shape == 'spiked_shape'
    
    len1 = A.dim(1)*unit;
    wid1 = A.dim(3)*unit;
    wid2 = A.dim(4)*unit;
    r1_x = A.r_x(1)*unit;
    rc1 = A.r(1)*unit;

    r2_x = A.r_x(2)*unit;
    rc2 = A.r(2)*unit;


    ECP1 = [q_x;a_y+wid1/2;q_x;a_y];
    ECP2 = [q_x-len1/2;a_y+wid1/2-wid2/2;r1_x+rc1;a_y];
    ECP3 = [q_x;a_y;q_x;a_y];
    ECP4 = [q_x;a_y;r1_x+rc1;a_y];
    ECP5 = [q_x-len1/2;a_y+wid1/2-wid2/2;r2_x+rc2;a_y];
    ECP6 = [q_x;a_y;r2_x+rc2;a_y];


    sig1 = 0;
    sig2 =0;
    sig3 = v_x;
    sig4 = 0;
    sig5 = 0;
    sig6 = 0;

    p_n1 = 0;
    p_n2 = 0;
    p_n3 = (m*g*h+F_elect*h);
    p_n4 = 0;
    p_n5 = 0;
    p_n6 = 0;

    p_t1 = 0;
    p_t2 = 0;
    p_t3 = mu*p_n3;
    p_t4 = 0;
    p_t5 = 0;
    p_t6 = 0;

    Con_wrench = [p_t1;p_t2;p_t3;p_t4;p_t5;p_t6];

    A.l(1:33,1) = -infty;
    A.l(34:75,1) = 0;
    A.u(1:75,1) = infty;
    La1 =[0;0;1;0;wid1/2]; %% assuming planar contact
    La2 =[0;0;0;0;0]; %% assuming planar contact
    La3 =[0;0.5;0.5;0;0]; %% assuming planar contact
    La4 =[0;0.5;0.5;0;0]; %% assuming planar contact
    La5 =[0;0;0;0;0]; %% assuming planar contact
    La6 =[0;0.5;0.5;0;0]; %% assuming planar contact

    A.fun = 'mcp_funjac_microrobot_spiked_2D_2cylinders_terrain';
    A.check = @mcp_funjac_microrobot_spiked_2D_2cylinders_terrain;

    A.Z = [nu;ECP1;ECP2;ECP3;ECP4;ECP5;ECP6;Con_wrench;sig1;sig2;sig3;sig4;sig5;sig6;La1;La2;La3;La4;La5;La6;p_n1;p_n2;p_n3;p_n4;p_n5;p_n6]; 

elseif A.shape == 'spiked_ended'
    len1 = A.dim(1)*unit;
    len2 = A.dim(2)*unit;
    wid1 = A.dim(3)*unit;
    wid2 = A.dim(4)*unit;
    
    r1_x = A.r_x(1)*unit;
    rc1 = A.r(1)*unit;

    r2_x = A.r_x(2)*unit;
    rc2 = A.r(2)*unit;
    
    ECP1 = [q_x-(len1/2-len2);wid2;r1_x+rc1;a_y];
    ECP2 = [q_x-(len1/2-len2);wid2;r2_x+rc2;a_y];
    
    ECP3 = [q_x+(len1/2-len2/2);a_y;q_x+(len1/2-len2/2);a_y];
    ECP4 = [q_x+(len1/2-len2/2);a_y;r1_x+rc1;a_y];
    ECP5 = [q_x+(len1/2-len2/2);a_y;r2_x+rc2;a_y];
    
    ECP6 = [q_x-(len1/2-len2/2);a_y;q_x-(len1/2-len2/2);a_y];
    ECP7 = [q_x-(len1/2-len2/2);a_y;r1_x+rc1;a_y];
    ECP8 = [q_x-(len1/2-len2/2);a_y;r2_x+rc2;a_y];
    
    sig1 = 0;
    sig2 = 0;
    sig3 = 0;
    sig4 = 0;
    sig5 = 0;
    sig6 = 0;
    sig7 = 0;
    sig8 = 0;
    
    p_n1 = 0;
    p_n2 = 0;
    p_n3 = 0;
    p_n4 = 0;
    p_n5 = 0;
    p_n6 = (m*g*h+F_elect*h);
    p_n7 = 0;
    p_n8 = 0;
    
    p_t1 = 0;
    p_t2 = 0;
    p_t3 = 0;
    p_t4 = 0;
    p_t5 = 0;
    p_t6 = mu*p_n6;
    p_t7 = 0;
    p_t8 = 0;
    
    Con_wrench = [p_t1;p_t2;p_t3;p_t4;p_t5;p_t6;p_t7;p_t8];
    
    A.l(1:43,1) = -infty;
    A.l(44:105,1) = 0;
    A.u(1:105,1) = infty;
    
    La1 =[0;0.7;0;0;0.5]; %% assuming planar contact
    La2 =[0;0.7;0;0;2]; %% assuming planar contact
    
    La3 =[0.5;0.5;0;0;0;0]; %% assuming planar contact
    La4 =[0;0.28;0;0;0;1.5]; %% assuming planar contact
    La5 =[0;0;0;0;0;2.92]; %% assuming planar contact
    
    La6 = [0.5;0.5;0;0;0;0]; %% assuming planar contact
    La7 = [0;0;0;0;0;0]; %% assuming planar contact
    La8 = [0;0;0;0;0;1.78]; %% assuming planar contact
    
    A.fun = 'mcp_funjac_microrobot_spiked_ended_2D_2cylinders_terrain';
    A.check = @mcp_funjac_microrobot_spiked_ended_2D_2cylinders_terrain;

    A.Z = [nu;ECP1;ECP2;ECP3;ECP4;ECP5;ECP6;ECP7;ECP8;Con_wrench;sig1;sig2;sig3;sig4;sig5;sig6;sig7;sig8;La1;La2;La3;La4;La5;La6;La7;La8;p_n1;p_n2;p_n3;p_n4;p_n5;p_n6;p_n7;p_n8]; 

    
end