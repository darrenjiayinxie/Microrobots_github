function [velocity] = velocity_results(A)
parameter =A.q(2,:);
N = size(parameter,2);
[scalar_max_displacement,index] = max(abs(A.q(2,:))); % maximum displacement from zero
vector_max_displacement = (A.q(2,index)/abs(A.q(2,index)))*scalar_max_displacement; % this code (along with the line above) accounts for the case
                                                                                    % where the robot might go up first before sliding down
velocity = (vector_max_displacement/(N*A.h))/1000; %[m/s]
% T = 1:N;
% plot(T*A.h,parameter);
% xlabel('Time (s)');
% ylabel('Horizontal Displacement (mm)');
% axis fill
end
