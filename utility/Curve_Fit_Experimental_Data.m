% Clean Up Previous Scripts/Values
clc
clear
close

% External Field Parameters
field_angle = 90;                                              % Angle of External Magnetic Field (°)

% Microrobot Properties
length = 0.8e-3;                                               % Length (m)
width = 0.4e-3;                                                % Width (m)
height = 0.1e-3;                                               % Height (m)
diagonal_length = sqrt(length^2+height^2);                     % Diagonal Length (m)
v = length*width*height;                                       % Volume (m^3)
m = 2750*(0.5e-3*width*height) + 1200*(0.3e-3*width*height);   % Mass (kg)
offset_angle = 35;                                             % Magnetic Alignment Offset Angle (°)

% Environment Parameters
g = 9.81;                                                      % Gravitational Constant (m/s^2)

% Experimental Data
B = 1e-3.*[3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10];                        % External Magnetic Field Strength (T) - Last Four Values Removed
theta = [9.819 18.662 24.039 32.754 38.494 43.365 44.444 46.042 47.657...
     48.731 51.085 50.425 51.892 50.721 51.87];                                  % Inclination Angle of Microrobot (°) - Last Four Values Removed
 
% Original Experimental Data
% B = 1e-3.*[3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 10.5 11 11.5 15];      % External Magnetic Field Strength (T)
% theta = [9.819 18.662 24.039 32.754 38.494 43.365 44.444 46.042 47.657...
%      48.731 51.085 50.425 51.892 50.721 51.87 51.564 53.341 53.384 55.231];    % Inclination Angle of Microrobot (°)

% Find Missing Parameters
fun =@(unknowns) norm(((diagonal_length/2)*(m*g + ...
    abs(unknowns(1)))*cosd(theta+atand(height/length)))...
     - B.*(v*abs(unknowns(2))*cosd(theta+offset_angle)));                        % Function to be Minimized (Based on Paper's Torque Equation)

unknowns = fminsearchbnd(fun,[5.48e-6,15000],[0, Inf],[0, Inf]);                 % Minimize the Function w/ Initial Guess and Boundary Conditions

force_elect = unknowns(1);                                                       % Electrostatic Force (N)
magnetization = unknowns(2);                                                     % Magnetization (A/m)

% Generate Model Data
theta_model = 0:0.01:90;
B_model = zeros(size(theta_model,2),1);
for i = 1:size(theta_model,2)
    B_model(i) = ((diagonal_length/2)*(m*g + abs(force_elect))*cosd(theta_model(i)...
        + atand(height/length))) / (v*abs(magnetization)*cosd(theta_model(i)+offset_angle));
end

% Plot Model Data Against Experimental Data
subplot(2,1,1)
plot(B,theta)
axis([0e-3 10e-3 0 90])
hold on
plot(B_model, theta_model)
title('Plot of Inclination Angle vs External Magnetic Field Strength')
xlabel('External Magnetic Field Strength (T)')
ylabel('Inclination Angle (°)')
legend('Experimental Data', 'Model Data')

% Find Error Between Model and Experimental Data
subplot(2,1,2)
error = zeros(size(theta,2),1);
for i = 1:size(theta,2)
    error(i) = ((diagonal_length/2)*(m*g + abs(force_elect))*cosd(theta(i)...
        + atand(height/length))) / (v*abs(magnetization)*cosd(theta(i)+offset_angle))...
        - B(i);
end
plot(B,error)
axis([0e-3 10e-3 -2e-3 Inf])
title('Plot of Error Between Model and Experiment Data')
xlabel('External Magnetic Field Strength (T)')
ylabel('Torque Error (N*m)')

% Print Results
fprintf('Electrostatic Force: %d [N]\n', force_elect);
fprintf('µTUM Magnetization: %d [A/m]\n', magnetization);