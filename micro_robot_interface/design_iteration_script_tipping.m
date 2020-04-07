% Overhead/Clean-up Functions
clear all
cd('..');
addpath('micro_robot_interface');
addpath('visualize');
addpath('utility_tipping');
addpath('visualize');
addpath('pathmexmaci64');
addpath('funjac/single_convex_contact_patches');
cd('micro_robot_interface')

% Simulation Parameters
n = 8;                      % number of iterations
freq = 0;                   % [Hz]
shape = 'spiked_shape';     % shape of robot

% Parameter Variation Data Set
len2 = [50 100 150 200 250 300 350 400 50 100 150 200 250 300 350 400];  %[microns]
wid1 = [50 100 150 200 250 300 350 400 400 350 300 250 200 150 100 50];  %[microns]


%% Tipping Data Collection

% Initialize Data Matrix and Iteration Variables
output = zeros(3,n);
n1 = 1;
n2 = 8;

for i = n1:n2
   % Print the Current Number of Iterations
   fprintf('Current Iteration: %i\n',i)
   
   % Reset Iteration Variables
   tilt_angle = 10;
   check = 1;
   
   % Record New Data Values
   output(1,i) = len2(i);
   output(2,i) = wid1(i);
   
   % Check When Robot Tips and Record Corresponding Angle
   while check == 1
        check = interface_micro_robot_tipping_function(shape,len2(i),wid1(i),tilt_angle,freq);
        if check == 1
            output(3,i) = tilt_angle;
            fprintf('    Angle: %.1f°\n',tilt_angle);
            tilt_angle = tilt_angle + 1;
        end
   end
end