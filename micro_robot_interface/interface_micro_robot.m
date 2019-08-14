
clear all
cd('..');
addpath('micro_robot_interface');
addpath('visualize');
addpath('utility');

addpath('funjac/single_convex_contact_patches');
addpath('pathmexmaci64');

cd('micro_robot_interface')
A = struct();

A.shape ='spiked_shape'; %'cuboid_shape', 'spiked_shape', 'spiked_ended', 'geckod_shape'
A = robot_inertia_parameters(A);


A.gravity = 9.8; %(m/s^2)
%% time step
A.h = 5e-4; %(s)

%% friction parameter
A.ellipsoid = [1 1 0.1]; % the choice of e_r (m) depends on the CM

A.material = 'paper';
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
A.M_r = 15000;% A/m old robot
%A.M_r = 51835;% A/m new robot
%A.phi = (27/180)*pi; % rad offset angle of old robot
A.phi = (0/180)*pi; % rad offset angle of new robot

A.B =20e-3; % T



% the frequency of the ratational magnetic field
A.fqn = 5; %Hz from 0Hz to 15 Hz


%% geometary of the inclined surface
A.theta = 0;
%A.theta = (45/180)*pi; % inclined angle


%% initial state

if A.shape == 'cuboid_shape'
     H = A.dim(3); 
elseif A.shape == 'spiked_shape'
     H = A.dim(3);
elseif A.shape == 'spiked_ended'
     H = A.dim(3)+2*A.dim(4);
elseif A.shape == 'geckod_shape'
     H = A.dim(3);
elseif A.shape == 'curved_shape'
     H = 2*(r1-D);
end
alpha = atan(H/A.dim(1));
theta_p = pi/2 - A.phi -alpha;

A.initial_q = [0;0;H/(2*cos(A.theta ));cos((A.theta +pi)/2);-sin((A.theta+pi)/2);0;0];
 
A.initial_v = [0;0;0;0;0;0];  %m/s

%% planner for the robot
A = planner(A);

A.n_i = 1;
A = NCP_micro_robot(A);
movie_microrobot(A);