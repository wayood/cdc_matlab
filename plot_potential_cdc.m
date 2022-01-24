load('potential_8_v1.mat');%軌道補正なし
load('potential_cdc_8_v1.mat');%軌道補正込み
load('potential_cruise_8.mat');%指示軌道ポテンシャル場
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
po_cruise_st=po_cruise;

%微分値とFFTによる前処理でポテンシャル遷移を導出
po_cruise=differential(po_cruise);
po_cdc=differential(po_cdc);
po=differential(po);
po_cruise=FFT(po_cruise);
po=FFT(po);
po_cdc=FFT(po_cdc);


%動的計画法による非線形マッチングの利用（軌道補正あり）
[cost,from]=DP_prepare(po_cruise,po_cdc);
[tmp_x,tmp_y]=DP_matching(cost,from);
[po_cruise,po_cdc]=GradientSR(po_cruise,po_cdc,tmp_x,tmp_y);

%補正後のLocal軌道
for i=1:length(po_cruise)
 sr_st_initial_cdc(i)=(po_cruise(i)-sum_po_cruise).^2;
 sr_st_cdc(i)=(po_cdc(i)-sum_po_cdc).^2;
 sr_st_cdc_con(i)=(po_cruise(i)-sum_po_cruise)*(po_cdc(i)-sum_po_cdc);
end

%動的計画法による非線形マッチングの利用（軌道補正なし）
[cost,from]=DP_prepare(po_cruise_st,po);
[tmp_x,tmp_y]=DP_matching(cost,from);
[po_cruise,po]=GradientSR(po_cruise_st,po,tmp_x,tmp_y);

for i=1:length(po_cruise)
 sr_st_initial_nocdc(i)=(po_cruise(i)-sum_po_cruise).^2;
 sr_st_nocdc(i)=(po(i)-sum_po).^2;
 sr_st_nocdc_con(i)=(po_cruise(i)-sum_po_cruise)*(po(i)-sum_po);
end

sr_cdc=sum(sr_st_cdc_con)/sqrt(sum(sr_st_initial_cdc))/sqrt(sum(sr_st_cdc));
sr=sum(sr_st_nocdc_con)/sqrt(sum(sr_st_initial_nocdc))/sqrt(sum(sr_st_nocdc));
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

function p_i=FFT(p_i)
     dt=0.01;%サンプル周期
     
     %サンプル数はとってきた配列数から2^nに近い値に設定し、リサンプリング
     n_i=fix(log2(length(p_i)));
     N_i=2.^n_i;
     p_i_sam=resample(p_i,N_i,length(p_i));
     fq_i=linspace(0,1.0/dt,N_i);
     
     
     %FFT開始
     p_i_fft=fft(p_i_sam);
    
     %高周波成分カット
     cutoff=fq_i(length(fq_i))*0.02;
     for i=1:length(fq_i)
         if fq_i(i)>cutoff
             p_i_fft(i)=0;
         end
     end
     %IFFTで元の信号に戻す
     p_i_ifft=real(ifft(p_i_fft))*2;
     p_i=resample(p_i_ifft,length(p_i),N_i);
end
%微分値を出力
function p_ini=differential(p_i)  
    for i=1:length(p_i)-1
        p_ini(i)=p_i(i+1)-p_i(i);
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

