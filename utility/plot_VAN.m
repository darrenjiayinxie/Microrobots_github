function plot_VAN(A)
VAN =A.VAN;
p_z = A.z(3,:)/(1e7);
N = size(VAN,2);
T = 1:N;
plot(T*A.h,VAN);
legend('Adhesive Force (N)')
xlabel('Time (s)');
ylabel('Force (N)')
end