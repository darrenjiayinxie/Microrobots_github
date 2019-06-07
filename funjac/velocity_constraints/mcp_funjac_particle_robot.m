function [F, J, domerr] = mcp_funjac_particle_robot(z, jacflag)

%% initialize
z = z(:);
F = [];
J = [];
domerr = 0;

%% obtain value of global variables
global h;

global q_old ;
q_xo = q_old(1);
q_yo = q_old(2);

global nu;
v_x =nu(1);
v_y =nu(2);

global n d;

%% unknown variables

q_x = z(1);
q_y = z(2);

c = z(3);



%% Newton_Euler equations

F(1) = q_x - q_xo - h*v_x-c*n(1);
F(2) = q_y - q_yo - h*v_y-c*n(2);

%% non-penetration 
F(3) = n(1)*q_x+n(2)*q_y+d;
%% Jacobian matrix (Using chain rules)

if (jacflag)
    J = [1,0,-n(1);
         0,1,-n(2);
         n(1),n(2),0];

    J = sparse(J);
end
end
