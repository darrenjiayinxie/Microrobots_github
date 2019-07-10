clear all
n = 20; % number of iterations
len2 = linspace(50,400,n); %[microns]
wid2 = linspace(50,500,n); %[microns]
output = zeros(1,n);

% CAN ONLY RUN ONE GATHERING SET UP AT ONE TIME (comment out the other)

% Speed data gathering
incline_angle = 0; % [degrees]
for i = 1:n
   i % print the current number of iterations
   output(i) = abs(simulate_micro_robot(len2(i), wid2(i), incline_angle));
end

% Incline climbing data gathering
% for i = 1:n
%    i % print the current number of iterations
%    incline_angle = 26;
%    check = 0;
%    while check >= 0
%         incline_angle % print the current incline angle
%         output(i) = check;
%         check = simulate_micro_robot(len2(i), wid2(i), incline_angle);
%         incline_angle = incline_angle + 0.125;
%    end
% end