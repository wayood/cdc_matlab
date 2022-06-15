%% 軌道補正アルゴリズム
global cdc_length;
cdc_length=0.5;%m単位で補正ポイントを設定

%% パラメータ設定
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
po_cdc=[];
slip=30;
glo_slip_x = 0;
glo_slip_y = 0;
dt=0.1;
drive_cdc=[];

for i=1:320
    ran_x=-50+100*rand;
    ran_y=100*rand;
    glo_rand_size(i)=0.3+1*rand;
    glo_obs(1,i)=ran_x;
    glo_obs(2,i)=ran_y;
end

[glo_gosa_obs]=gosamodel(glo_obs,p);%経路生成時のLM座標
graph(glo_gosa_obs,p,glo_rand_size);
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
for io=1:length(glo_obs(1,:))
  en_plot_red(glo_obs(:,io).',glo_rand_size(io));
  hold on;
end
for j=1:length(wp(1,:))
  plot(wp(1,j),wp(2,j),'g:o','MarkerSize',10);
  hold on;
end

glo_start=p.start;
%%  逐次補正プログラム
%軌道補正経路 リアルタイム

%初期化
i=1;
count=1;
sr=0;
p_cdc=0;
p_cd=0;
p_init=0;
p_i=0;
ang_wp = pi/2;
v=5;

%ナビゲーション
while 1
   p.start=move_rover(p.start,wp(:,i));
   l=len(wp(:,i).',p.start.');
   if count>1
     delete(b);%逐次的に処理を削除(現在のLM座標(障害物))
     delete(w);
   end
   if count==1
    [me_gosa_obs]=gosa_hozon(glo_gosa_obs);
   end
   [h,~]=size(me_gosa_obs);
   if  h == 3
      me_gosa_obs(3,:)=[];
   end
   if l<1.0
       %wpでのsaferateを計算→前原さん、宮本さん参照
       p_init(i)=potential(glo_obs,wp_init(:,i),glo_rand_size);
       p_i=sum(p_init)/i;%ポテンシャル場の平均走行前軌道
       p_cdc(i)=potential(glo_obs,wp(:,i),glo_rand_size);
       p_cd=sum(p_cdc)/i;%走行後のパテンシャル場の平均
       if i == length(wp(1,:))
            disp("Finish");
            break;
       end
       i=i+1;
   end
   [cur_obs]=gosa_move(me_gosa_obs,p.start.',ang_wp,v);
   [ang_wp,sen_num,up_obs]=sensor_range(cur_obs,p.start,ang_wp);
   [rand_size,gosa_obs]=sensor_judge(glo_gosa_obs,sen_num,glo_rand_size);
   [b,w]=animation(up_obs,wp,p,i,rand_size);
   %ここでアニメーションが完成  
   [wp,k,mat_er,plan_er]=correction(up_obs,gosa_obs,p,i);
   count=count+1;
   pause(0.001);
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
   plot(two_point(1,:),two_point(2,:),'-b');
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
%% 誤差モデル計算　
%距離による関係性を考えた。
%これはDEMのステレオデータによる誤差を主に考えた。
function [up_obs]=gosamodel(obs,p)
 for i=1:length(obs(1,:))
    r(i)=len(obs(:,i).',p.start.');
    if r(i)<1
        r(i)=0;
    end
    l(i)=18.5*10^-2*r(i)^2*0.01;
    x_ran=-0.01+0.02*rand; %キャリブレーションや分解能での誤差を考える
    y_ran=-0.01+0.02*rand;
    up_obs(1,i)=obs(1,i)+x_ran;
    up_obs(2,i)=obs(2,i)+l(i)+y_ran;
end
    up_obs(3,:)=1;
end

%誤差をロボット進行方向に考えて表示
function [up_obs]=gosa_move(obs,start,ang_wp,v)
 global glo_slip_x;
 global glo_slip_y;
 sliprate(ang_wp,v);
 for i=1:length(obs(1,:))
    r(i)=len(obs(:,i).',start.');
    if r(i)<1
        r(i)=0;
    end
    l(i)=18.5*10^-2*r(i)^2*0.01;
    error=l(i);
    error_x=cos(ang_wp)*error;
    error_y=sin(ang_wp)*error;
    up_obs(1,i)=error_x+obs(1,i)+glo_slip_x;
    up_obs(2,i)=error_y+obs(2,i)+glo_slip_y;
    
 end
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
function [awp,k,mat_er,plan_er]=correction(lm_current,lm_first,p,i)
    global wp_init;
    A=lm_current*pinv(lm_first);
    wp(3,:)=1;
    A(3,1)=0;
    A(3,2)=0;
    [h,~]=size(lm_current);
    if h == 2
        lm_current(3,:) = 1;
    end
    global A_n;
    if isempty(A_n)==0
     [mat_er,plan_er] = A_matrix(A,lm_current,lm_first,A_n);
    else
        mat_er=0;
        plan_er=0;
    end
    A_n=A;
   for j=i:length(wp_init(1,:))
       wp(1,j)=wp_init(1,j);
       wp(2,j)=wp_init(2,j);
   end
   for i=1:length(wp(1,:))
      awp(:,i)=A()*wp(:,i);
   end
   awp(1,:)=awp(1,:);
   awp(2,:)=awp(2,:);
   awp(3,:)=[];
   k=cond(A,2);
end

%% A行列解析
%カリタニさん参照
%ネットに条件数とかの同じような説明ありそちらも参照可能
function [mat_er,plan_er]=A_matrix(A,lm_c,lm_f,A_n)
  mat_er = cond(A)*norm(A*lm_f-lm_c)/norm(A*lm_f);
  plan_er = cond(A_n)*norm(A_n-A)/norm(A_n);
end

%% センサ検知前位置のグラフを表示
function graph(glo_gosa_obs,p,size)
  for i=1:length(glo_gosa_obs(1,:))
    en_plot_blue(glo_gosa_obs(:,i).',size(i));
  end
    plot(p.start(1,1),p.start(2,1),'b:.','MarkerSize',5);
    hold on;
    grid on;
    xlabel('x[m]')
    ylabel('y[m]')
    xlim([-50 50]);
    ylim([0 100]);
end

%% ロボットが観測せず移動している間の軌跡のプロット
function fx(af_start,bef_start)
    A=[bef_start(1,1) 1;
        af_start(1,1) 1];
    B=[bef_start(2,1);
        af_start(2,1)];
    X=linsolve(A,B);
    a=X(1,1);
    b=X(2,1);
    if af_start(2,1)>bef_start(2,1)
      for y=bef_start(2,1):0.01:af_start(2,1)
         plot((y-b)/a,y,'b:.','MarkerSize',3);
         hold on;
       end
    else
      for y=af_start(2,1):0.01:bef_start(2,1)
         plot((y-b)/a,y,'b:.','MarkerSize',3);
         hold on;
       end
    end
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
   range_base=50;

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
 length_x=goal(1,1)-start(1,1);
 length_r=sqrt((goal(1,1)-start(1,1))^2+(goal(2,1)-start(2,1))^2);
 ang=acos(length_x/length_r);
end

%% 二点間のプロット
function f_twopoint(af_start,bef_start)
  A=[bef_start(1,1) 1;
     af_start(1,1) 1];
  B=[bef_start(2,1);
     af_start(2,1)];
  X=linsolve(A,B);
  a=X(1,1);
  b=X(2,1);
  y=linspace(bef_start(2,1),af_start(2,1));
  plot((y-b)/a,y,'r:.','MarkerSize',3);
  hold on;
end

%% アニメーション　円の作成
%下三つは円の塗りつぶし
function en_plot_blue(glo_obs,size)
 [x,y]=circle(glo_obs(1,1),glo_obs(1,2),size);
 fill(x,y,'b');
 hold on;
end

function a=en_plot_red(glo_obs,size)
 [x,y]=circle(glo_obs(1,1),glo_obs(1,2),size);
 a=fill(x,y,'r');
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

