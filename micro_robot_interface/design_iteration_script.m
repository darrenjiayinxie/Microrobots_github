clear all
cd('..');
addpath('micro_robot_interface');
addpath('visualize');
addpath('pathmexmaci64');
addpath('funjac/single_convex_contact_patches');
addpath('utility');

cd('micro_robot_interface')
n = 20; % number of iterations
len2 = linspace(50,400,n); %[microns]
wid2 = linspace(50,500,n); %[microns]

% CAN ONLY RUN ONE GATHERING SET UP AT ONE TIME (comment out the other)
%output = zeros(1,n);

% Speed data gathering
% incline_angle = 0; % [degrees]
% for i = 1:n
%    i % print the current number of iterations
%    output(i) = abs(simulate_micro_robot(len2(i), wid2(i), incline_angle));
% end

%Incline climbing data gathering
output = zeros(2,n);
n1 = 1;
n2 = 20;
for i = n1:n2
   i % print the current number of iterations
   incline_angle = 26;
   check = 0;
   out = [0;0];
   while check >= 0
        incline_angle; % print the current incline angle
       
        output(:,i) = out;
        out = simulate_micro_robot(len2(i), wid2(i), incline_angle);
        incline_angle = incline_angle + 0.125;
        check = out(1);
   end
end