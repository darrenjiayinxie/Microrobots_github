function check = interface_micro_robot_tipping_function(shape,len2,wid1,tilt_angle,freq)

wid1 = wid1*2 + 150;                                    % correct definition of wid1
check = 1;                                              % initialize check variable

A = struct();

A.shape = shape;
if A.shape == 'cuboid_shape'
    A.dim=[800e-6 100e-6 600e-6];                       % (m) length width height
    A.density = 2.1688e3;                               % kg/m^3
    
elseif A.shape == 'spiked_shape'
    A.dim=[800e-6 len2*1e-6 wid1*1e-6 150e-6 400e-6];   % (m) len1 len2 wid1 wid2 heg
    A.density = 2.1688e3;                               % kg/m^3
    
elseif A.shape == 'spiked_ended'
    A.dim=[800e-6 len2*1e-6 150e-6 wid1*1e-6 400e-6];   % (m) len1 len2 wid2 wid1 heg
    A.density = 2.1688e3;                               % kg/m^3

end

A = robot_inertia_parameters(A);


A.gravity = 9.8; % (m/s^2)
%% time step
A.h = 5e-4; % (s)

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
A.M_r = 15000;% A/m old robot
%A.M_r = 51835;% A/m new robot
%A.phi = (27/180)*pi; % rad offset angle of old robot
A.phi = (0/180)*pi; % rad offset angle of new robot

A.B =20e-3; % T



% the frequency of the ratational magnetic field
A.fqn = freq; %Hz from 0Hz to 15 Hz


%% geometary of the inclined surface
%A.theta = 0;
A.theta = (tilt_angle/180)*pi; % inclined angle


%% initial state

if A.shape == 'cuboid_shape'
     H = A.dim(3); 
elseif A.shape == 'spiked_shape'
     H = A.dim(3);
elseif A.shape == 'spiked_ended'
     H = A.dim(3)+2*A.dim(4);
end
A.initial_q = [0;0;H/2;cos(pi/2);-sin(pi/2);0;0];
 
A.initial_v = [0;0;0;0;0;0];  %m/s => [v_x,v_y,v_z,w_x,w_y,w_z]

%% units
% unit for meter
A.unit = 1e3; % 1 - m, 10 - dm, 100 - cm, 1000 - mm ,1e6 um 
% unit for mass
A.unit_mass = 1e3; % 1-kg, 1e3 - g;%% planner for the robot
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

eul = A.q(4:7,:);
eul = transpose(eul);
eul = quatern2euler(eul);

% movie_microrobot_tipping(A);

for i = 1:size(eul,1)
    if abs(eul(i,2)) > 0.5
        check = 0;
    elseif tilt_angle > 60
        check = 0;
    end
end

function euler = quatern2euler(q)
%QUATERN2EULER Converts a quaternion orientation to ZYX Euler angles
%
%   q = quatern2euler(q)
%
%   Converts a quaternion orientation to ZYX Euler angles where phi is a
%   rotation around X, theta around Y and psi around Z.
%
%   For more information see:
%   http://x-io.co.uk/quaternions/
%
%	Date          Author          Notes
%	27/09/2011    SOH Madgwick    Initial release

    R(1,1,:) = 2.*q(:,1).^2-1+2.*q(:,2).^2;
    R(2,1,:) = 2.*(q(:,2).*q(:,3)-q(:,1).*q(:,4));
    R(3,1,:) = 2.*(q(:,2).*q(:,4)+q(:,1).*q(:,3));
    R(3,2,:) = 2.*(q(:,3).*q(:,4)-q(:,1).*q(:,2));
    R(3,3,:) = 2.*q(:,1).^2-1+2.*q(:,4).^2;

    phi = atan2(R(3,2,:), R(3,3,:) );
    theta = -atan(R(3,1,:) ./ sqrt(1-R(3,1,:).^2) );
    psi = atan2(R(2,1,:), R(1,1,:) );

    euler = [phi(1,:)' theta(1,:)' psi(1,:)'];
end

end