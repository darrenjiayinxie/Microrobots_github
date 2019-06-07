function [X,Y,Z] = mashgrid_cuboid(l,w,h,H,q_b)

d_w = w/1;
d_l = l/1;
d_h = h/1;
% bottom face in body frame
[X1,Y1] = meshgrid(-w/2:d_w:w/2,-l/2:d_l:l/2);
Z1 = X1.*0 + Y1.*0 -h/2;


% top face in body frame
[X2,Y2] = meshgrid(-w/2:d_w:w/2,-l/2:d_l:l/2);
Z2 = X2.*0 + Y2.*0 +h/2;

% left face in body frame
[X3,Z3] = meshgrid(-w/2:d_w:w/2,-h/2:d_h:h/2);
Y3 = X3.*0 + Z3.*0 -l/2;


% right face in body frame
[X4,Z4] = meshgrid(-w/2:d_w:w/2,-h/2:d_h:h/2);
Y4 = X4.*0 + Z4.*0 +l/2;


% front face in body frame
[Y5,Z5] = meshgrid(-l/2:d_l:l/2,-h/2:d_h:h/2);
X5 = Y5.*0 + Z5.*0 +w/2;

% back face in body frame
[Y6,Z6] = meshgrid(-l/2:d_l:l/2,-h/2:d_h:h/2);
X6 = Y6.*0 + Z6.*0 -w/2;


X = [X1;X2;X3;X4;X5;X6];
Y = [Y1;Y2;Y3;Y4;Y5;Y6];
Z = [Z1;Z2;Z3;Z4;Z5;Z6];

H_b = [1,0,0,q_b(1);
       0,1,0,q_b(2);
       0,0,1,q_b(3);
       0,0,0,1];


[m,n] = size(X);
Q_x = zeros(m,n);
Q_y = zeros(m,n);
Q_z = zeros(m,n);

for i = 1:m
   
    A = H*H_b*[X(i,:);Y(i,:);Z(i,:);ones(1,n)];
    %A(1:3,:) = R*A(1:3,:);
    Q_x(i,:) = A(1,:);
    Q_y(i,:) = A(2,:);
    Q_z(i,:) = A(3,:);
end
X = Q_x;
Y = Q_y;
Z = Q_z;

end