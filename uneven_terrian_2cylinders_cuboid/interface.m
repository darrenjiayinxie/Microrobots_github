%% notation
clear all;
close all;
clc


addpath('data');
addpath('visualize');
cd('..');
addpath('funjac/single_convex_contact_patches');
addpath('pathmexmaci64');

cd('uneven_terrian_2cylinders_cuboid')

A = struct();

% approximate dimension of the object 
A.shape ='cuboid_shape';
A.cylinder =1;
A.dim=[800e-6 400e-6 100e-6]; %(m) length width height
A.density = 2.1688e3; %kg/m^3
A.mass = 1.6071e-7; %(kg)
A.V_m = 2.9e-11; % m^3d


A.gravity = 9.8; %(m/s^2)
% time step
A.h = 5e-4; %(s)

% friction model
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



%A.M_r = 15000;% A/m old robot
A.M_r = 51835;% A/m new robot
%A.phi = (27/180)*pi; % rad offset angle of old robot
A.phi = (0/180)*pi; % rad offset angle of new robot

A.B =10e-3; % T
%A.B =0;





% the frequency of the ratational magnetic field
A.fqn = 5; %Hz from 0Hz to 15 Hz





% geometary of the inclined surface
%A.theta = 0;
A.theta = (0/180)*pi; % inclined angle

% geometry of the terrain
delta_1 = 0.05e-3; %distance between robot and first bump
A.rc_1 =  0.3e-3;
A.r1_y =- (delta_1 + A.rc_1 +A.dim(1)/2 + A.dim(3));
A.r1_z = 0;

delta_2 = 0.5e-3; %distance between first bump and second bump
A.rc_2 =  0.6e-3;
A.r2_y = - (delta_2 + A.rc_2 +A.rc_1) +A.r1_y;
A.r2_z = 0;


% initial state
H = A.dim(3); 
   
A.initial_q = [0;0;H/(2*cos(A.theta ));cos((A.theta +pi)/2);-sin((A.theta+pi)/2);0;0];

  
A.initial_v = [0;0;0;0;0;0];  %m/s

% unit for meter
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