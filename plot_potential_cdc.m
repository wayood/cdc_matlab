load('potential.mat');%軌道補正なし
load('potential_cdc.mat');%軌道補正込み
Potential_Field(glo_obs,glo_rand_size);
hold on;
%% 三次元的にポテンシャル場で評価
for i=1:length(drive(1,:))
    z(i)=potential(glo_obs,drive(:,i),glo_rand_size);
    plot3(drive(1,i)+11,drive(2,i),z(i),'o','Color','b');
    hold on;
end
for j=1:length(drive_cdc(1,:))
    z_cdc(j)=potential(glo_obs,drive_cdc(:,j),glo_rand_size);
    plot3(drive_cdc(1,j)+11,drive_cdc(2,j),z_cdc(j),'o','Color','r');
    hold on;
end

%% 補正によるSafeRateを表現
p_init=sum(z)/i;
p_cdc=sum(z_cdc)/j;

for i=1:length(p_init)
 sr_st1(i)=(z(i)-p_init).^2;
 sr_st2(i)=(z_cdc(i)-p_cdc).^2;
 sr_st(i)=(z(i)-p_init)*(z_cdc(i)-p_cdc);
end
sr=sum(sr_st)/sqrt(sum(sr_st1))/sqrt(sum(sr_st2));
T=table(sr);
T=T(1,1);
fig = uifigure;
uit = uitable(fig,'Data',T);
hold on;

%% ポテンシャル場を二次元プロット
pause();
hold off;
i=1:length(drive(1,:));
plot(i,z(i));
hold on;
i=1:length(drive_cdc(1,:));
plot(i,z_cdc(i));
hold on;
grid on;
xlabel('x[m]');
ylabel('potential');



function po=potential(obs,move,size)
    po=0;
    for i=1:obs(1,:)+1
     l=norm(obs(:,i).'-move.');
     if l < size(i)
       p=(3-l.^2/size(i).^2)/2;
     else
       p=size(i)/l;
     end
       po=po+p;
    end
end