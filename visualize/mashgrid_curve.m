function [X,Y,Z] = mashgrid_curve(r1,r2,D,wid,alpha,H,q_b)


[X1,Y1] = meshgrid(-wid/2:wid/20:wid/2,-sin(alpha)*r1:2*sin(alpha)*r1/20:sin(alpha)*r1);
Z1 = -sqrt(r1^2-Y1.^2)+D;

[X2,Y2] = meshgrid(-wid/2:wid/20:wid/2,-sin(alpha)*r2:2*sin(alpha)*r2/20:sin(alpha)*r2);
Z2 = -sqrt(r2^2-Y2.^2)+D;

[X3,Y3] = meshgrid(-wid/2:wid/20:wid/2,sin(alpha)*r2:sin(alpha)*(r1-r2)/20:sin(alpha)*r1);
Z3 = D-(cos(alpha)/sin(alpha))*Y3;

[X4,Y4] = meshgrid(-wid/2:wid/20:wid/2,-sin(alpha)*r1:sin(alpha)*(r1-r2)/20:-sin(alpha)*r2);
Z4 = D+(cos(alpha)/sin(alpha))*Y4;


X = [X1;X4;X2;X3];
Y = [Y1;Y4;Y2;Y3];
Z = [Z1;Z4;Z2;Z3];

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