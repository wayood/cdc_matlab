%% 軌道補正アルゴリズム
clear all;
global cdc_length;
cdc_length=0.5;%m単位で補正ポイントを設定

%% 初期宣言
p.start=[0;0];
global glo_obs;
global glo_gosa_obs;
global glo_rand_size;
global drive_cdc;
global dt;
global slip;
global glo_slip_x;
global glo_slip_y;
global po_cdc;
global lm_cur_1;
global Path_analysis;
Path_analysis =[];
lm_cur_1 = [];
glo_slip_x = 0;
glo_slip_y = 0;
drive_cdc=[];
ang_wp = pi/2;
range_base=20;
i=1;

%% 障害物
while 1
    ran_x=-20+40*rand;
    ran_y=100*rand;
    [ang,l] = cart2pol(ran_x,ran_y);
    glo_rand_size(i)=0.3+0.5*rand;
    glo_obs(1,i)=ran_x;
    glo_obs(2,i)=ran_y;   
    
    if l <= range_base && (7*pi)/36 <= ang && ang <= (29*pi)/36
        glo_rand_size(i)=0.3+0.5*rand;
        glo_obs(1,i)=ran_x;
        glo_obs(2,i)=ran_y;
        if i == 50
            break;
        end
        i=i+1;
    end    
end

[glo_gosa_obs]=gosamodel(glo_obs,p);%経路生成時のLM座標
gosa_plot = graph(glo_gosa_obs,p,glo_rand_size);

%wpを手動で設定
[x,y]=ginput;
wp=[x.';
    y.'];
global wp_init;
wp_init = wp;

f_twopoint(wp(:,1),p.start());
for i=1:length(wp(1,:))-1
  f_twopoint(wp(:,i+1),wp(:,i));
end

hold on;

for j=1:length(wp(1,:))
  plot(wp(1,j),wp(2,j),'g:o','MarkerSize',10);
  hold on;
end

glo_start=p.start;
%%  逐次補正プログラム
%軌道補正経路 リアルタイム

%初期化
i=1;
po_cdc=[];
slip=30;
dt=0.1;
count=1;
sr=0;
p_cdc=0;
p_cd=0;
p_init=0;
p_i=0;
v=5;

obs = ob_round(glo_gosa_obs,glo_rand_size);
%ナビゲーション
while 1
   ang_wp = angular(wp(:,i),p.start);
   start = [p.start;
            ang_wp];
   p.start=move_rover(p.start,wp(:,i));
   l=len(wp(:,i).',p.start.');
   if count>1
     delete(b);%逐次的に処理を削除(現在のLM座標(障害物))
     delete(w);
     delete(L);
   end
   if count==1
    [me_gosa_obs]=gosa_hozon(glo_gosa_obs);
   end
   [h,~]=size(me_gosa_obs);
   if  h == 3
      me_gosa_obs(3,:)=[];
   end
   if l<0.5
       %wpでのsaferateを計算→前原さん、宮本さん参照
       p_init(i)=potential(glo_obs,wp_init(:,i),glo_rand_size);
       p_i=sum(p_init)/i;%ポテンシャル場の平均走行前軌道
       p_cdc(i)=potential(glo_obs,wp(:,i),glo_rand_size);
       p_cd=sum(p_cdc)/i;%走行後のパテンシャル場の平均
       if i == length(wp(1,:))
            disp("Finish");
            [x,y]=ginput;
            wp_add=[x.';
                    y.'];
             pl_wp=[wp(:,end) wp_add];
             plot(pl_wp(1,:),pl_wp(2,:),'-r','LineWidth',2);
             plot(pl_wp(1,:),pl_wp(2,:),'g:o','MarkerSize',10);
             wp = [wp wp_add];
             wp_init = wp;
             glo_gosa_obs(1,:) = glo_gosa_obs(1,:) +glo_slip_x;
             glo_gosa_obs(2,:) = glo_gosa_obs(2,:) +glo_slip_y;
             delete(gosa_plot);
             glo_gosa_obs(3,:) = [];
             for plt_cou = 1:length(glo_gosa_obs(1,:))
                [x,y]=circle(glo_gosa_obs(1,plt_cou),glo_gosa_obs(2,plt_cou),glo_rand_size(plt_cou));
                gosa_plot(plt_cou)=fill(x,y,'b');
                hold on;
             end
             glo_gosa_obs(3,:) = 1;
             start = [p.start;
                      ang_wp];
            
       end
       i=i+1;
   end
   [cur_obs]=gosa_move(me_gosa_obs,p.start.',ang_wp,v);
   [ang_wp,sen_num,up_obs]=sensor_range(cur_obs,p.start,ang_wp);
   [rand_size,gosa_obs]=sensor_judge(glo_gosa_obs,sen_num,glo_rand_size);
   [b,w]=animation(up_obs,wp,p,i,rand_size);
   %ここでアニメーションが完成  
   [wp,k,mat_er,plan_er]=correction(up_obs,gosa_obs,start,i);
   L = lm_line(up_obs,gosa_obs);
   count=count+1;
   pause(0.1);
end
%% 逐次的に進むローバの角度算出と走行距離の決定
function start=move_rover(start,wp)
   global cdc_length;
   r=cdc_length;%進む距離を調整
   ang=angular(wp,start);
   bef_start=start;
   start(1,1)=r*cos(ang)+start(1,1);
   start(2,1)=r*sin(ang)+start(2,1);
   two_point=[start bef_start];
   plot(two_point(1,:),two_point(2,:),'-b','LineWidth',2);
end
%% ポテンシャル場で評価
function po=potential(obs,move,size)
    po=0;
    for i=1:length(obs(1,:))
     l=len(obs(:,i).',move.');
     if l < size(i)
       p=(3-l.^2/size(i).^2)/2;
     else
       p=size(i)/l;
     end
     po=po+p;
    end
end

function kill = lm_line(LM_current,LM_first)
    LM_first(3,:) = [];
    for i = 1:length(LM_first(1,:))
        LM = [LM_current(:,i) LM_first(:,i)];
        kill(i)=plot(LM(1,:),LM(2,:),'-g','LineWidth',1.5);
        hold on;
    end
end
%% 誤差モデル計算　
%距離による関係性を考えた。
%これはDEMのステレオデータによる誤差を主に考えた。
function [up_obs]=gosamodel(obs,p)
 for i=1:length(obs(1,:))
    [ang,leng] = cart2pol(obs(1,i),obs(2,i));
    r(i)=len(obs(:,i).',p.start.');
    if r(i)<1
        r(i)=0;
    end
    x_ran=-0.01+0.02*rand; %キャリブレーションや分解能での誤差を考える
    y_ran=-0.01+0.02*rand;
    if leng >= 20
        l=18.5*10^-2*20^2*0.01;
        up_obs(1,i)=obs(1,i)+x_ran;
        up_obs(2,i)=obs(2,i)+l+y_ran;
        continue;
    end
    l=18.5*10^-2*r(i)^2*0.01;
    up_obs(1,i)=obs(1,i)+x_ran;
    up_obs(2,i)=obs(2,i)+l+y_ran;
end
    up_obs(3,:)=1;
end

%誤差をロボット進行方向に考えて表示
function [up_obs]=gosa_move(obs,start,ang_wp,v)
 global glo_slip_x;
 global glo_slip_y;
 sliprate(ang_wp,v);
up_obs(1,:)=obs(1,:)+glo_slip_x;
up_obs(2,:)=obs(2,:)+glo_slip_y;
up_obs(3,:)=1;
end

%% 逐次的なアニメーション
%逐次的に軌道補正を表示
function [b,w]=animation(up_obs,wp,p,i,size)
  for j=i:length(wp(1,:))
  w(j)=plot(wp(1,j),wp(2,j),'g:o','MarkerSize',10);
  end
  %plot(p.start(1,1),p.start(2,1),'r:.','MarkerSize',3);
  %hold on;
  for j=1:length(up_obs(1,:))
  b(j)=en_plot_orange(up_obs(:,j).',size(j));
  end
end

%距離を計算
function l=len(a,b)
l=norm(a-b);
end

%% 軌道補正計算
function [awp,k,mat_er,plan_er]=correction(lm_current,lm_first,start,i)
    global wp_init;
    global lm_cur_1;
    A=lm_current*pinv(lm_first);
    wp_init(3,:)=1;
    A(3,1)=0;
    A(3,2)=0;
    [h,~]=size(lm_current);
    if h == 2
        lm_current(3,:) = 1;
    end
    global A_n;
    if isempty(A_n)==0
     [mat_er,plan_er] = A_matrix(A,lm_current,lm_first,A_n,lm_cur_1);
    else
        mat_er=0;
        plan_er=0;
    end
    lm_cur_1 = lm_current;
    A_n=A;
    awp=A*wp_init;
    awp(3,:)=[];
    wp_init(3,:) = [];
    k=cond(A,2);
end

%% A行列解析
%カリタニさん参照

function [mat_er,plan_er]=A_matrix(A,LM_current,LM_first,A_n,LM_t_1)
    global Path_analysis;
    mat_er = cond(A)*norm(A*LM_first-LM_current)/norm(A*LM_first);
    plan_er = cond(A_n)*norm(A_n-A)/norm(A_n);
    [h,~] = size(LM_t_1);
    [~,Z_first,V_first] = svd(LM_first);
    [~,Z_current,V_current] = svd(LM_current);
    VTRate_spatial = sum(dot(V_current,V_first))/3;
    CNRate_spatial = cond(Z_current)/cond(Z_first); 
    if h == 3
        [~,Z_t_1,~] = svd(LM_t_1);
        %VTRate_seque = sum(dot(V_current,V_t_1))/3;
        CNRate_seque = cond(Z_current)/cond(Z_t_1);
        Path_analysis_vir = [mat_er,plan_er,VTRate_spatial,CNRate_spatial,CNRate_seque];
        Path_analysis = [Path_analysis;Path_analysis_vir];
        file = sprintf("A_matrix.mat");
        save(file,"Path_analysis");
        fprintf("A matrix disperation --> cond %f,%f\n",CNRate_spatial,VTRate_spatial);
    end
end

%% センサ検知前位置のグラフを表示
function [b] = graph(glo_gosa_obs,p,size)

  for i=1:length(glo_gosa_obs(1,:))
    b(i) = en_plot_blue(glo_gosa_obs(:,i).',size(i));
  end
    plot(p.start(1,1),p.start(2,1),'b:.','MarkerSize',5);
    hold on;
    grid on;
    xlabel('x[m]')
    ylabel('y[m]')
    xlim([-15 15]);
    ylim([0 30]);
end


%% 経路生成時のLM座標の視野角判定
function [rand_size,cur_gosa_obs]=sensor_judge(gosa_obs,sen_num,glo_rand_size)
   for i=1:length(gosa_obs(1,:))
      if sen_num(i)==1
          gosa_obs(:,i)=[-1;-1;-1];
          glo_rand_size(i)=0;
      end
   end
    idx = gosa_obs(1,:)==-1 & gosa_obs(2,:) == -1 & gosa_obs(3,:)==-1;
    idy = glo_rand_size(1,:)==0;
    rand_size=glo_rand_size(~idy);
    cur_gosa_obs=gosa_obs(:,~idx);
end

%% 視野角を考慮
function [ang_wp,sen_num,sen_obs]=sensor_range(obs,start,ang)
   %視野の射程距離
   range_base=20;

   if obs(3,:)==1
       obs(3,:)=[];
   end
   
   ang_90=pi/2-ang;
   b=start(2,1)-start(1,1)*tan(ang_90);
   for i=1:length(obs(1,:))
     [ang_wp,range_max,range_min,range_s,range_ln]=siya(obs(:,i),start,ang);
    
     if range_s>range_min && range_s<range_max && range_ln<range_base
        sen_num(i)=0;
     else
        obs(:,i)=[-1;-1];
        sen_num(i)=1;
     end
   end
      idx = obs(1,:)== -1 & obs(2,:) == -1;
      sen_obs = obs(:,~idx);
end

%% 視野角計算
function   [ang,range_wpbase_max,range_wpbase_min,range_s,range_l]=siya(obs,start,ang)
   range_wpbase_min=ang-(11*pi/36);
   range_wpbase_max=ang+(11*pi/36);
   range_l=sqrt((obs(1,1)-start(1,1))^2+(obs(2,1)-start(2,1))^2);
   range_x=obs(1,1)-start(1,1);
   range_y=obs(2,1)-start(2,1);
   range_s=atan2(range_y,range_x);
end

function ang=angular(goal,start)
 [ang,~] = cart2pol(goal(1,1)-start(1,1),goal(2,1)-start(2,1));
end

%% 二点間のプロット
function f_twopoint(af_start,bef_start)
    start = [af_start bef_start];
    plot(start(1,:),start(2,:),'-r','LineWidth',2);
    hold on;
end

%% アニメーション　円の作成
%円の塗りつぶし
function b=en_plot_blue(glo_obs,size)
 [x,y]=circle(glo_obs(1,1),glo_obs(1,2),size);
 b=fill(x,y,'b');
 hold on;
end

function b=en_plot_orange(glo_obs,size)
 orange=[0.9500 0.6250 0];%オレンジ
 [x,y]=circle(glo_obs(1,1),glo_obs(1,2),size);
 b=fill(x,y,orange);
 hold on;
end

%円の関数
function [r_x,r_y]=circle(x,y,r)
 t=linspace(0,2*pi,100);
 r_x=r*cos(t)+x;
 r_y=r*sin(t)+y;
end

%蓄積誤差を保存
function [up_obs]=gosa_hozon(obs)
 global glo_slip_x;
 global glo_slip_y;
 if obs(3,:) == 1
     obs(3,:)=[];
 end
 for i=1:length(obs(1,:))
    up_obs(1,i)=obs(1,i)+glo_slip_x;
    up_obs(2,i)=obs(2,i)+glo_slip_y;
 end
 up_obs(3,:)=1;
end
%% 障害物の円を座標格納
function obs=ob_round(cur_obs,r)
for i=1:length(r)
    [x,y]=circle(cur_obs(1,i),cur_obs(2,i),r(i));
    if i==1
        obs=[x;y];
    else
    obs=[obs(1,:) x;obs(2,:) y];
    end
end
end
%% スリップ率導入
function sliprate(ang,v)
 global dt; 
 global glo_slip_x;
 global glo_slip_y;
 global slip;
 s=rand*slip/100;
 v_real=(1-s)*v;
 slip_length=(v_real-v)*dt;
 x=slip_length*cos(ang);
 y=slip_length*sin(ang);
 glo_slip_x =  glo_slip_x + x;
 glo_slip_y =  glo_slip_y + y;
end

