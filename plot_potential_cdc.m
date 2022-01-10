load('potential_8_v1.mat');%軌道補正なし
load('potential_cdc_8_v1.mat');%軌道補正込み
x=100;
y=100;
Potential_Field(glo_obs,glo_rand_size,x,y);
hold on;
%% 三次元的にポテンシャル場で評価
%ここでは真の障害物に対しての安定性を補正ありなしで判別

for i=1:length(drive(1,:))
    z(i)=potential(glo_obs,drive(:,i),glo_rand_size);
    plot3(drive(1,i)+51,drive(2,i),z(i),'o','Color','b');
    hold on;
end
s1=sum(z)/i;

for j=1:length(drive_cdc(1,:))
    z_cdc(j)=potential(glo_obs,drive_cdc(:,j),glo_rand_size);
    plot3(drive_cdc(1,j)+51,drive_cdc(2,j),z_cdc(j),'o','Color','r');
    hold on;
end
s2=sum(z_cdc)/j;
pause();

%% 補正によるSafeRateを表現
%SafeRateをGlobalpathでのポテンシャル遷移を計算
glo_gosa_obs(3,:)=[];
for i=1:length(path(:,1))
    z(i)=potential(glo_gosa_obs,path(i,:).',glo_rand_size);
end
p_init=sum(z)/i;

%補正後のLocal軌道を
if length(po_cdc) <= length(z)
    num=length(po_cdc);
else
    num=length(z);
end

for i=1:length(z)
 sr_st_initial(i)=(z(i)-p_init).^2;
end

for i=1:num
 sr_st_cdc(i)=(po_cdc(i)-sum_po_cdc).^2;
 sr_st_cdc_con(i)=(z(i)-p_init)*(po_cdc(i)-sum_po_cdc);
end

if length(po) <= length(z)
    num=length(po);
else
    num=length(z);
end

for i=1:num
 sr_st_nocdc(i)=(po(i)-sum_po).^2;
 sr_st_nocdc_con(i)=(z(i)-p_init)*(po(i)-sum_po);
end
sr_cdc=sum(sr_st_cdc_con)/sqrt(sum(sr_st_initial))/sqrt(sum(sr_st_cdc));
sr=sum(sr_st_nocdc_con)/sqrt(sum(sr_st_initial))/sqrt(sum(sr_st_nocdc));
SR=[sr_cdc,sr];
T=array2table(SR,'VariableNames',{'SR_cdc','SR'});
fig = uifigure;
uit = uitable(fig,'Data',T);
hold on;

%% ポテンシャル場を二次元プロット
pause();
hold off;
i=1:length(drive(1,:));
plot(i,z(i),'b');
hold on;
i=1:length(drive_cdc(1,:));
plot(i,z_cdc(i),'r');
hold on;
grid on;
xlabel('x[m]');
ylabel('potential');
save("potential_evaluation","sr","sr_cdc","s1","s2");


function po=potential(obs,move,size)
    po=0;
    for i=1:length(obs(1,:))
     l=norm(obs(:,i).'-move.');
     if l < size(i)
       p=(3-l.^2/size(i).^2)/2;
     else
       p=size(i)/l;
     end
       po=po+p;
    end
end