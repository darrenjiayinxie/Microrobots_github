 
clear all
cd('..');
addpath('micro_robot_interface');
addpath('visualize');
addpath('utility');

addpath('funjac/single_convex_contact_patches');
addpath('pathmexmaci64');

cd('micro_robot_interface')
A = struct();

A.shape ='spiked_shape';
if A.shape == 'cuboid_shape'
    A.dim=[800e-6 400e-6 100e-6]; %(m) length width height
    A.density = 2.1688e3; %kg/m^3
    
elseif A.shape == 'spiked_shape'
    A.dim=[400e-6 125e-6 325e-6 75e-6 200e-6]; %(m) len1 len2 wid1 wid2 heg
    %A.density = 2.1688e3; %su-8 %kg/m^3
    A.density = 2.6684e3;  % PDMS
    A.err = [pi*4/180 pi*4/180]; %derr_1 derr_2
elseif A.shape == 'spiked_ended'
    A.dim=[800e-6 100e-6 150e-6 125e-6 400e-6]; %(m) len1 len2 wid1 wid2 heg
    A.density = 2.1688e3; %kg/m^3

end

A = robot_inertia_parameters(A);


A.gravity = 9.8; %(m/s^2)
%% time step
A.h = 5e-4; %(s)

%% friction parameter
A.ellipsoid = [1 1 0.1]; % the choice of e_r (m) depends on the CM

A.material = 'PDMSS';
if A.material == 'paper'
    A.cof = 0.3; % coefficient of friction
    A.van_constant = 3.7148; %N/m^2
    A.F_elect = 3.2022e-6; %N
elseif A.material == 'alumi'
    A.cof = 0.5359;
    A.van_constant = 26.1771; %N/m^2
    A.F_elect = 0;
elseif A.material == 'PDMSS'
    A.cof = 0.8;
    A.van_constant = 340.115; %N/m^2
    A.F_elect = 2.2022e-7;
end
%


%% magnetic parameters
A.M_r = 18661;% A/m old robot
%A.M_r = 51835;% A/m new robot
%A.phi = (27/180)*pi; % rad offset angle of old robot
A.phi = (20/180)*pi; % rad offset angle of new robot

A.B =20e-3; % T



% the frequency of the ratational magnetic field
A.fqn = 1; %Hz from 0Hz to 15 Hz


%% geometary of the inclined surface
%A.theta = 0;
A.theta = (60/180)*pi; % inclined angle


%% initial state

if A.shape == 'cuboid_shape'
     H = A.dim(3); 
elseif A.shape == 'spiked_shape'
     H = A.dim(3);
elseif A.shape == 'spiked_ended'
     H = A.dim(3)+2*A.dim(4);
end
A.initial_q = [0;0;H/(2*cos(A.theta ));cos((A.theta +pi)/2);-sin((A.theta+pi)/2);0;0];
 
A.initial_v = [0;0;0;0;0;0];  %m/s => [v_x,v_y,v_z,w_x,w_y,w_z]

%% units
% unit for meter
A.unit = 1e3; % 1 - m, 10 - dm, 100 - cm, 1000 - mm ,1e6 um 
% unit for mass
A.unit_mass = 1e3; % 1-kg, 1e3 - g;%% planner for the robot
A = planner(A);

A.n_i = 1;
A = NCP_micro_robot(A);
movie_microrobot(A);