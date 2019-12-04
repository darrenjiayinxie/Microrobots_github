function [F1, F2] = compute_F_elect(A)
cof = A.cof;
m = A.mass;
g = A.gravity;
theta1 = (45/180)*pi;
theta2 = (70/180)*pi;

F1 = m*g*sin(theta1)/cof - m*g*cos(theta1);
F2 = m*g*sin(theta2)/cof - m*g*cos(theta2);
end