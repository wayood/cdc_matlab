load('potential_9_v1.mat');%軌道補正なし
load('potential_cdc_9_v1.mat');%軌道補正込み
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


[z,po_cdc]=differential(z,po_cdc);
hold off;
i=1:length(z(1,:));
plot(i,z(i),'b');
hold on;
i=1:length(po_cdc(1,:));
plot(i,po_cdc(i),'r');
ylim([-0.02 0.02]);
xlim([0 length(z)]);
hold on;
grid on;
pause();
[z,po_cdc]=FFT(z,po_cdc);
i=1:length(z(1,:));
plot(i,z(i),'r');
hold on;
i=1:length(po_cdc(1,:));
plot(i,po_cdc(i),'b');
hold on;
ylim([-0.02 0.02]);
xlim([0 length(z)]);
grid on;
pause();
[cost,from]=DP_prepare(z,po_cdc);
[tmp_x,tmp_y]=DP_matching(cost,from);
[z,po_cdc]=GradientSR(z,po_cdc,tmp_x,tmp_y);
pause();
hold off;
i=1:length(z(1,:));
plot(i,z(i),'b');
hold on;
i=1:length(po_cdc(1,:));
plot(i,po_cdc(i),'r');
hold on;
grid on;
p_init=sum(z)/i;

%補正後のLocal軌道を

for i=1:length(z)
 sr_st_initial(i)=(z(i)-p_init).^2;
 sr_st_cdc(i)=(po_cdc(i)-sum_po_cdc).^2;
 sr_st_cdc_con(i)=(z(i)-p_init)*(po_cdc(i)-sum_po_cdc);
end

[cost,from]=DP_prepare(z,po);
[tmp_x,tmp_y]=DP_matching(cost,from);
[z,po]=GradientSR(z,po,tmp_x,tmp_y);

for i=1:length(z)
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

%% 関数系

%ポテンシャル場の生成
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

function [z,po_cdc]=FFT(z,po_cdc)
     cutoff=1;
     z=fft(z);
     po_cdc=fft(po_cdc);
     for i=1:length(z)
         if z(i)>cutoff
             z(i)=0;
         end
     end
     for i=1:length(po_cdc)
         if po_cdc(i)>cutoff
             po_cdc(i)=0;
         end
     end
     z=ifft(z);
     po_cdc=ifft(po_cdc);
end
%微分値を出力
function [p_ini,p_cdc]=differential(p_i,p_t)
    
    for i=1:length(p_i)-1
        p_ini(i)=p_i(i+1)-p_i(i);
    end

    for i=1:length(p_t)-1
        p_cdc(i)=p_t(i+1)-p_t(i);
    end
end

%最小コストとマップの作成
function [cost,from]=DP_prepare(p_i,p_t)
cost=zeros(length(p_i),length(p_t));
from=zeros(length(p_i),length(p_t));
for i=1:length(p_i)
    for j=1:length(p_t)
        d(i,j)=norm(p_i(i)-p_t(j));
        if i==1 && j==1
            cost(i,j)=d(i,j);
        end
        if i>1
            cost(i,1)=cost(i-1,1)+d(1,j);
            from(i,1)=1;
        end
        if j>1
            cost(1,j)=cost(1,j-1)+d(i,1);
            from(1,j)=2;
        end
        if i>1 && j>1
            A=[cost(i-1,j)+d(i,j) cost(i-1,j-1)+2*d(i,j) cost(i,j-1)+d(i,j)];
            [M,I]=min(A);
            cost(i,j)=M;
            if I==1
                from(i,j)=0;
            elseif I==2
                from(i,j)=1;
            else 
                from(i,j)=2;
            end
        end
    end
end

end

%最短経路を探索
function [tmp_x,tmp_y]=DP_matching(cost,from)
    i=length(cost(:,1));
    j=length(cost(1,:));
    leng=i+j;
    for k=1:leng
        if from(i,j)==0
            i=i-1;
        elseif from(i,j)==1
            i=i-1;
            j=j-1;
        else
            j=j-1;
        end
        if i<1
            i=1;
        end
        if j<1
            j=1;
        end
        if i==1 && j==1
            tmp_x(k)=i;
            tmp_y(k)=j;
            break;
        end 
        tmp_x(k)=i;
        tmp_y(k)=j;
    end
end

%対応するポテンシャル場対して割り振り
function [z,po_cdc]=GradientSR(p_i,p_t,tmp_i,tmp_t)
     i=length(tmp_i):-1:1;
     j=1:length(tmp_i);
     tmp_i(j)=tmp_i(i);
     i=length(tmp_i):-1:1;
     j=1:length(tmp_i);
     tmp_t(j)=tmp_t(i);
     for i=1:length(tmp_i)
         z(i)=p_i(tmp_i(i));
         po_cdc(i)=p_t(tmp_t(i));
     end
end

