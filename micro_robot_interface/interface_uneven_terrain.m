%% notation
clear all;
close all;
clc

cd('..');
addpath('funjac/single_convex_contact_patches');
addpath('pathmexmaci64');
addpath('visualize');
addpath('utility_uneven_terrain');
cd('micro_robot_interface');

A = struct();

%% approximate dimension of the object 
A.shape ='cuboid_shape';

if A.shape == 'cuboid_shape'
    A.dim=[800e-6 400e-6 100e-6]; %(m) length width height
    A.density = 2.1688e3; %kg/m^3
    
%elseif A.shape == 'spiked_shape'
 %   A.dim=[800e-6 200e-6 650e-6 150e-6 300e-6]; %(m) len1 len2 wid1 wid2 heg
 %   A.density = 2.1688e3; %kg/m^3
    
%elseif A.shape == 'spiked_ended'
 %  A.dim=[800e-6 100e-6 150e-6 125e-6 400e-6]; %(m) len1 len2 wid1 wid2 heg
 %  A.density = 2.1688e3; %kg/m^3

end
A.cylinder =1;
A = robot_inertia_parameters(A);

A.gravity = 9.8; %(m/s^2)
%% time step
A.h = 5e-4; %(s)

%% friction model
A.ellipsoid = [1 1 0.1]; % the choice of e_r (m) depends on the CM

A.material = 'alumi';
if A.material == 'paper'
    A.cof = 0.3; % coefficient of friction
    A.van_constant = 3.7148; %N/m^2
    A.F_elect = 3.2022e-6; %N
elseif A.material == 'alumi'
    A.cof = 0.5359;
    A.van_constant = 26.1771; %N/m^2
    A.F_elect = 0;
end
%


%% magnetic parameters
A.M_r = 51835;% A/m new robot
A.phi = (0/180)*pi; % rad offset angle of new robot
A.B =10e-3; % T
% the frequency of the ratational magnetic field
A.fqn = 10; %Hz from 0Hz to 15 Hz





%% geometary of the inclined surface
%A.theta = 0;
A.theta = (0/180)*pi; % inclined angle


%% initial state

if A.shape == 'cuboid_shape'
     H = A.dim(3); 
%elseif A.shape == 'spiked_shape'
 %    H = A.dim(3);
%elseif A.shape == 'spiked_ended'
 %    H = A.dim(3)+2*A.dim(4);
end
A.initial_q = [0;0;H/(2*cos(A.theta ));cos((A.theta +pi)/2);-sin((A.theta+pi)/2);0;0];
 
A.initial_v = [0;0;0;0;0;0];  %m/s


%% geometry of the terrain
N = 5; % number of bumps


[r_y,r_z,r] = generate_obstacles(A,N);
A.r_y = r_y;
A.r_z = r_z;
A.r = r';


%% unit for meter
A.unit = 1e3; % 1 - m, 10 - dm, 100 - cm, 1000 - mm ,1e6 um 
% unit for mass
A.unit_mass = 1e3; % 1-kg, 1e3 - g;


% planner for the input 
A = planner(A);

A.n_i = 1;
success =0;
while success ==0
    try
        A = NCP_micro_robot(A);
        success = 1;
    catch
        success = 0;
    end
end
    

%v = compute_velocity(A);
movie_microrobot_uneven_terrain(A);