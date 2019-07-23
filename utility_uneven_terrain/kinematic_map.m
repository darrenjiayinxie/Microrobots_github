function [q,nu] = kinematic_map(q_old,z,h)

%% updating the position
q = zeros(7,1);

nu=z(1:6);

q(1:3) = q_old(1:3)+h*nu(1:3);

%% updating the orientation

q0 = q_old(4);
q1 = q_old(5);
q2 = q_old(6);
q3 = q_old(7);

E = [-q1 q0 -q3 q2;
 -q2 q3 q0 -q1;
 -q3 -q2 q1 q0];

q_term = q_old(4:7)+h*0.5*E'*nu(4:6);
q_norm = (q_term(1)^2+q_term(2)^2+q_term(3)^2+q_term(4)^2)^(1/2); % normalize the Quaternion

q(4) = q_term(1)/q_norm;
q(5) = q_term(2)/q_norm;
q(6) = q_term(3)/q_norm;
q(7) = q_term(4)/q_norm;

end