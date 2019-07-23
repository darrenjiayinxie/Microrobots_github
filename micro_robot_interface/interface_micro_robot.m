
clear all
cd('..');
addpath('micro_robot_interface');
addpath('visualize');


addpath('funjac/single_convex_contact_patches');
addpath('pathmexmaci64');

cd('micro_robot_interface')
A = struct();

%% approximate dimension of the object 
A.shape ='spiked_shape';
if A.shape == 'cuboid_shape'
    A.dim=[800e-6 400e-6 100e-6]; %(m) length width height
    A.density = 2.1688e3; %kg/m^3
    A.mass = 1.6071e-7; %(kg)
    A.V_m = 2.9e-11; % m^3
% 
%     A.density = 2.1688e3; %kg/m^3
%     A.mass = 6.94e-8; %(kg)
%     A.V_m = 3.2e-11; % m^3
elseif A.shape == 'spiked_shape'
    A.dim=[800e-6 200e-6 650e-6 150e-6 300e-6]; %(m) len1 len2 wid1 wid2 heg
    A.density = 2.1688e3; %kg/m^3
%     A.mass = 1.1061e-7; %(kg)
%     A.V_m = 3.6e-11; % m^3
    A.mass = 6.94e-8; %(kg)
    A.V_m = 3.2e-11; % m^3
elseif A.shape == 'spiked_ended'
    A.dim=[800e-6 100e-6 150e-6 125e-6 300e-6]; %(m) len1 len2 wid1 wid2 heg
    A.density = 2.1688e3; %kg/m^3
%     A.mass = 9.4343e-8; %(kg)
%     A.V_m = 3.6e-11; % m^3
    A.mass = 6.94e-8; %(kg)
    A.V_m = 3.2e-11; % m^3
elseif A.shape == 'geckod_shape'
    A.dim=[800e-6 400e-6 270e-6 150e-6]; %(m) length width height
    A.density = 2.1688e3; %kg/m^3
    A.mass = 9.4343e-8; %(kg)
    A.V_m = 3.6e-11; % m^3
elseif A.shape == 'curved_shape'
    A.dim=[1250e-6 1100e-6 400e-6 (25/180)*pi]; %(m) r_1 r_2 wid theta
    r1 = A.dim(1);
    r2 = A.dim(2);
    wid = A.dim(3);
    alpha = A.dim(4);
    D = ((r1^3-r2^3)/(r1^2-r2^2))*(2*sin(alpha)/(3*alpha));
    A.density = 2.1688e3; %kg/m^3
%     A.V_m = (pi*r1^2*wid-pi*r2^2*wid)*(alpha*180/pi)/180; % m^3
%     A.mass = A.V_m*A.density; %(kg)
    A.mass = 6.94e-8; %(kg)
    A.V_m = 3.2e-11; % m^3
end




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

%% unit for meter
A.unit = 1e3; % 1 - m, 10 - dm, 100 - cm, 1000 - mm ,1e6 um 
% unit for mass
A.unit_mass = 1e3; % 1-kg, 1e3 - g;


%% planner for the robot
A = planner(A);

A.n_i = 1;
A = NCP_micro_robot(A);
movie_microrobot(A);