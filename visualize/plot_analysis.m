
T = 1:2000;
load('cuboid.mat');
z_cuboid = A.z;
q_cuboid = A.q;
load('spiked_end.mat');
z_spi_end = A.z;
q_spi_end = A.q;
load('spiked_shape.mat');
z_spi_shape = A.z;
q_spi_shape = A.q;
plot(T*A.h,q_cuboid(2,:),T*A.h,q_spi_end(2,:),T*A.h,q_spi_shape(2,:));

legend('cuboid','spiked end','spiked shape');