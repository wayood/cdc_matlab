load('A_matrix.mat');
t=1:1:numel(k_A);
subplot(2,1,1); plot(t,Matrix_error(t))
ylim([0 1]);
ylabel('Mapping Error Rate');
xlabel('cdc point(整数)');
grid on;
title('Mapping Error Rate')

t=1:1:numel(k_A);
subplot(2,1,2); plot(t,k_A(t))
ylim([1 100]);
ylabel('condition number');
xlabel('cdc point(整数)');
grid on;
title('condition number')

for i=1:length(p_init)
 sr_st1(i)=(p_init(i)-p_i).^2;
 sr_st2(i)=(p_cdc(i)-p_cd).^2;
 sr_st(i)=(p_init(i)-p_i)*(p_cdc(i)-p_cd);
end
sr=sum(sr_st)/sqrt(sum(sr_st1))/sqrt(sum(sr_st2));

T=table(sr);
T=T(1,1);
fig = uifigure;
uit = uitable(fig,'Data',T);
hold on;