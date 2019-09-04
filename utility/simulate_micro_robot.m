function [output] = simulate_micro_robot(len2,wid2,incline_angle)

%%  dimensions and density of the robot of the robot
A = struct();
A.shape ='spiked_shape';
if A.shape == 'cuboid_shape'
    A.dim=[800e-6 400e-6 100e-6]; %(m) length width height
    A.density = 2.1688e3; %kg/m^3
 
elseif A.shape == 'spiked_shape'
    wid1 = wid2*2+150;
    A.dim=[800e-6 len2*1e-6  wid1*1e-6 150*1e-6  300e-6]; %(m) len1 len2 wid1 wid2 heg
    A.density = 2.1688e3; %kg/m^3

elseif A.shape == 'spiked_ended'
    A.dim=[800e-6 len2*1e-6 150e-6 wid2*1e-6 300e-6]; %(m) len1 len2 wid1 wid2 heg
    A.density = 2.1688e3; %kg/m^3

end



%% calculate the volume and mass and moment of inertia of the robot 
A = robot_inertia_parameters(A);


A.gravity = 9.8; %(m/s^2)
%% time step
A.h = 5e-4; %(s)

%% friction parameter
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
%A.M_r = 15000;% A/m old robot
A.M_r = 51835;% A/m new robot
%A.phi = (27/180)*pi; % rad offset angle of old robot
A.phi = (0/180)*pi; % rad offset angle of new robot

A.B =20e-3; % T



% the frequency of the rotational magnetic field
A.fqn = 10; %Hz from 0Hz to 15 Hz


%% geometry of the inclined surface
A.theta = (incline_angle/180)*pi; % [radians]
%A.theta = (45/180)*pi; % inclined angle


%% initial state

if A.shape == 'cuboid_shape'
     H = A.dim(3); 
elseif A.shape == 'spiked_shape'
     H = A.dim(3);
elseif A.shape == 'spiked_ended'
     H = A.dim(3)+2*A.dim(4);
end


A.initial_q = [0;0;H/(2*cos(A.theta ));cos((A.theta +pi)/2);-sin((A.theta+pi)/2);0;0];
 
A.initial_v = [0;0;0;0;0;0];  %m/s

%% unit for meter
A.unit = 1e3; % 1 - m, 10 - dm, 100 - cm, 1000 - mm ,1e6 um 
% unit for mass
A.unit_mass = 1e3; % 1-kg, 1e3 - g;


%% planner for the robot
A = planner(A);

A.n_i = 1;
success = 0;
while success == 0
    try
        A = NCP_micro_robot(A); % if path solver couldn't find a solution, continue anyways and output '-1'
        success = 1;
    catch
        success = 0;
        disp('failed, try again');
    end
end
%v = compute_velocity(A);
%movie_microrobot(A);


if incline_angle == 0
    max_velocity = velocity_results(A); % output the maximum velocity if inputted incline angle is 0 degrees
    output = max_velocity;
else
    max_velocity = velocity_results(A);
    if max_velocity >= 0
        output = [-1;max_velocity];                    % output '-1' if the robot went backwards and slipped off the incline
    else
        output = [incline_angle;max_velocity];         % output the incline angle if the robot succeeded in climbing it
    end
end
