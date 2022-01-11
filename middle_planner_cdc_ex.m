%% 読み込みとプロット
clear all;
global glo_obs;
global glo_gosa_obs;
global glo_rand_size;
global drive_cdc;
global dt;
global slip;
global glo_slip_x;
global glo_slip_y;
global po_i;
po_i=1;
slip=0;
glo_slip_x = 0;
glo_slip_y = 0;
dt=0;
load("path_interp_8.mat");
path=drive_cdc;
drive_cdc=[];
p.start=[0;0];


graph(glo_gosa_obs,p,glo_rand_size);

for io=1:length(glo_obs(1,:))
  en_plot_red(glo_obs(:,io).',glo_rand_size(io));
  hold on;
end

if length(path)>=1
    plot(path(1,:),path(2,:),'-r','MarkerSize',3);
    hold on;
end
hold on;

%% 曲率

[x,y]=ginput;
wp=[x.';
    y.'];
%{
wp_x=path(1,:).';
wp_y=path(2,:).';
p_wp=[wp_x,wp_y];
[L2,R2,K2] = curvature(p_wp);
quiver(wp_x,wp_y,K2(:,1),K2(:,2));
hold on;
j=1;

con_wp_x=[];
con_wp_y=[];
for i=1:length(wp_x)
    if abs(K2(i,1)) > 1 || abs(K2(i,2)) > 1 || i==length(wp_x)
     con_wp_x(j)=wp_x(i);
     con_wp_y(j)=wp_y(i);
     j=j+1;
    end
end
   
wp = [con_wp_x;con_wp_y];
%}
save("wp_8_v2.mat","wp");
for i=1:length(wp(1,:))
  kill(i)=plot(wp(1,i),wp(2,i),'g:o','MarkerSize',10);
  hold on;
end

point=twopoint(wp(:,1),p.start);
kill_point=zeros(100,5000);
for i=1:length(wp(1,:))-1
   kill_point(:,i)=twopoint(wp(:,i+1),wp(:,i));
end
kill_point(length(wp(:,1)),:)=[];
    
pause();

global wp_init;
wp_init = wp;
global p_i;
p_i=1;
%%  逐次補正プログラム
%軌道補正経路 リアルタイム

%初期化
i=1;
count=1;
obs=ob_round(glo_gosa_obs,glo_rand_size);

%ナビゲーション
while 1
   [p.start,up_obs,gosa_obs]=DynamicWindowApproachSample_k(p.start.',wp(:,i).',obs.',path.');
   l=len(wp(:,i).',p.start.');
   if count==1
   delete(kill);
   delete(point);
   delete(kill_point);
   end
   if count>1
     %逐次的に処理を削除(現在のLM座標(障害物))
     delete(kill_p);
     delete(w);
   end
   i=i+1;
   %[b,w]=animation(up_obs,wp,p,i,rand_size);
   %ここでアニメーションが完成
   if len(wp(:,length(wp(1,:))).',p.start.')<0.7 
       fx(wp(:,i),p.start);
       disp("Finish");
       break;
   end
   [wp,k,mat_er,plan_er]=correction(up_obs,gosa_obs,p,i-1);
   k_A(count)=k;
   Matrix_error(count)=mat_er;
   Plan_error(count)=plan_er;
   count=count+1;
   c_path=[];
   path=[];
   for j=i-1:length(wp(1,:))-1
       kill_p_ex=twopoint(wp(:,j+1),wp(:,j));
       c_path=path_correct(wp(:,j+1),wp(:,j));
       path=[path c_path];
       w(j)=plot(wp(1,j),wp(2,j),'g:o','MarkerSize',10);
       kill_p(j-i+2,:)=kill_p_ex;
   end
   pause(0.001);
end


%% ポテンシャル場で評価
function po=potential(obs,move,size)
    po=0;
    if obs(3,:)==1
        obs(3,:)=[];
    end
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
 s=randi(60)/100;
 v_real=(1-s)*v;
 slip_length=(v_real-v)*dt;
 x=slip_length*cos(ang);
 y=slip_length*sin(ang);
 glo_slip_x =  glo_slip_x + x;
 glo_slip_y =  glo_slip_y + y;
end

%% 逐次的なアニメーション
%逐次的に軌道補正を表示
function [b,w]=animation(up_obs,wp,p,i,size)
  for j=i:length(wp(1,:))
  w(j)=plot(wp(1,j),wp(2,j),'g:o','MarkerSize',10);
  end
  plot(p.start(1,1),p.start(2,1),'r:.','MarkerSize',3);
  hold on;
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

%% 誤差を明確化ラインプロット
function kill=cdc_obs_line(lm_first,lm_current)
for i=1:length(lm_first(1,:))
    kill(i)=twopoint(lm_current(:,i),lm_first(:,i));
end
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
function [ang_wp,sen_num,sen_obs]=sensor_range(obs,start,wp)
   %視野の射程距離
   range_base=50;
   if obs(3,:)==1
       obs(3,:)=[];
   end
   for i=1:length(obs(1,:))
     [ang_wp,range_wpbase_max,range_wpbase_min,range_s,range_ln]=siya(obs(:,i),start,wp);
     if start(2,1)>obs(2,i) 
        obs(:,i)=[-1;-1];
        sen_num(i)=1;%視野に入るかの判定
      elseif range_s>range_wpbase_min && range_s<range_wpbase_max && range_ln<range_base
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
function   [ang_wp,range_wpbase_max,range_wpbase_min,range_s,range_l]=siya(obs,start,wp)
   range_l=sqrt((obs(1,1)-start(1,1))^2+(obs(2,1)-start(2,1))^2);
   range_x=obs(1,1)-start(1,1);
   range_l1=sqrt((wp(1,1)-start(1,1))^2+(wp(2,1)-start(2,1))^2);
   range_x1=wp(1,1)-start(1,1);
   range_wpbase_min=acos(range_x1/range_l1)-(11*pi/36);
   range_wpbase_max=acos(range_x1/range_l1)+(11*pi/36);
   ang_wp=acos(range_x1/range_l1);
   range_s=acos(range_x/range_l);
end

function ang=angular(goal,start)
 length_x=goal(1,1)-start(1,1);
 length_r=sqrt((goal(1,1)-start(1,1))^2+(goal(2,1)-start(2,1))^2);
 ang=acos(length_x/length_r);
end

%% 二点間のプロット
function t=twopoint(af_start,bef_start)

  a=(bef_start(2,1)-af_start(2,1))/(bef_start(1,1)-af_start(1,1));
  b=bef_start(2,1)-a*bef_start(1,1);
  x=linspace(bef_start(1,1),af_start(1,1));
  y=linspace(bef_start(2,1),af_start(2,1));
  
  if bef_start(1,1)-af_start(1,1)==0
   b=bef_start(1,1);
   t=plot(b,y,'g:.','MarkerSize',3);
   hold on;
  elseif bef_start(2,1)-af_start(2,1)==0
   t=plot(x,b,'g:.','MarkerSize',3);
   hold on;
  else
   t=plot((y-b)/a,y,'g:.','MarkerSize',3);
   hold on;
  end 
  
end

%% 軌道の保存
function t=path_correct(af_start,bef_start)

  a=(bef_start(2,1)-af_start(2,1))/(bef_start(1,1)-af_start(1,1));
  b=bef_start(2,1)-a*bef_start(1,1);
  x=linspace(bef_start(1,1),af_start(1,1));
  y=linspace(bef_start(2,1),af_start(2,1));
  
  if bef_start(1,1)-af_start(1,1)==0
   b=bef_start(1,1);
   t=[b;y];
  elseif bef_start(2,1)-af_start(2,1)==0
   t=[x;b];
  else
   t=[(y-b)/a;y];
  end 
  
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

%% local plan
function [s,up_obs,gosa_obs] = DynamicWindowApproachSample_k(start,goal,obstacle,path)

x=[start pi/2 0 0]';%ロボットの初期状態[x(m),y(m),yaw(Rad),v(m/s),ω(rad/s)]
global glo_obs;
global glo_gosa_obs;
global glo_rand_size;
global drive_cdc;
obstacleR=0.5;%衝突判定用の障害物の半径
global dt; 
global po_i;
dt=0.1;%刻み時間[s]

%ロボットの力学モデル
%[最高速度[m/s],最高回頭速度[rad/s],最高加減速度[m/ss],最高加減回頭速度[rad/ss],
% 速度解像度[m/s],回頭速度解像度[rad/s]]
Kinematic=[1.0,toRadian(20.0),0.2,toRadian(50.0),0.01,toRadian(1)];

%評価関数のパラメータ [heading,dist,velocity,predictDT,path]
evalParam=[0.1,0.2,0.1,3.0];

%シミュレーション結果
result.x=[];
%tic;
% Main loop

for i=1:5000
%DWAによる入力値の計算
if i==1
    [u,traj]=DynamicWindowApproach(x,Kinematic,goal,evalParam,obstacle,obstacleR,path);
else
    [u,traj]=DynamicWindowApproach(x,Kinematic,goal,evalParam,obs,obstacleR,path);
end

x=f(x,u);%運動モデルによる移動

%シミュレーション結果の保存
result.x=[result.x; x'];
start=[x(1),x(2)];
s_x=[x(1);x(2)];
drive_cdc=[drive_cdc s_x];

if i==1
    [me_gosa_obs]=gosa_hozon(glo_gosa_obs);
end
[ang_wp,sen_num,cur_obs]=sensor_range(me_gosa_obs,start.',goal.');
[up_obs]=gosa_move(cur_obs,start.',x(3),u(1,1));
[rand_size,gosa_obs]=sensor_judge(glo_gosa_obs,sen_num,glo_rand_size);
ob=ob_round(up_obs,rand_size);
obs=ob.';
%kill_line=cdc_obs_line(gosa_obs,up_obs);
po_cdc(po_i)=potential(up_obs,s_x,rand_size);
sum_po_cdc=sum(po_cdc)/po_i;
po_i=po_i+1;
save('potential_cdc_8_v2.mat','glo_obs','glo_gosa_obs','glo_rand_size','drive_cdc','po_cdc','sum_po_cdc');
if i>1
    delete(d_q);
    delete(d_g);
    delete(d_tr);
    delete(b);
end

for j=1:length(up_obs(1,:))
  b(j)=en_plot_orange(up_obs(:,j).',rand_size(j));
end


%ゴール判定
if norm(x(1:2)-goal')<3.0
    disp('Arrive Goal!!');
    s=[x(1);x(2)];
    delete(b);
    break;
end

if  i > 30 && abs(result.x(length(result.x(:,1)),1) - result.x(length(result.x(:,1))-30,1)) < 1.0 && abs(result.x(length(result.x(:,1)),2) - result.x(length(result.x(:,1))-30,2)) < 1.0 
    obstacleR=0.1;
end

if  i > 40 && abs(result.x(length(result.x(:,1)),1) - result.x(length(result.x(:,1))-40,1)) < 1.0 && abs(result.x(length(result.x(:,1)),2) - result.x(length(result.x(:,1))-40,2)) < 1.0 
    evalParam(1)=2.0;    
    evalParam(2)=0.2;
    evalParam(5)=0.3;
end

if  i > 100 && abs(result.x(length(result.x(:,1)),1) - result.x(length(result.x(:,1))-100,1)) < 1.0 && abs(result.x(length(result.x(:,1)),2) - result.x(length(result.x(:,1))-100,2)) < 1.0 
    disp('Skip Waypoint');
    s=[x(1);x(2)];
    delete(b);
    break;
end
if i>1    
    delete(d_x);   
end
    
    %====Animation====
 ArrowLength=0.5;%矢印の長さ
 %ロボット
 d_q=quiver(x(1),x(2),ArrowLength*cos(x(3)),ArrowLength*sin(x(3)),'ok');
 hold on;
 d_x=plot(result.x(:,1),result.x(:,2),'-b');
 hold on;
 d_g=plot(goal(1),goal(2),'*r');
 hold on;
 %探索軌跡表示
 if ~isempty(traj)
    for it=1:length(traj(:,1))/5
       ind=1+(it-1)*5;
       d_tr(it)=plot(traj(ind,:),traj(ind+1,:),'-g');
       hold on;
    end
 end
 drawnow;
end
%toc
end
%movie2avi(mov,'movie.avi');
 

function [u,trajDB]=DynamicWindowApproach(x,model,goal,evalParam,ob,R,path)
%DWAによる入力値の計算をする関数

%Dynamic Window[vmin,vmax,ωmin,ωmax]の作成
Vr=CalcDynamicWindow(x,model);
%評価関数の計算
[evalDB,trajDB]=Evaluation(x,Vr,goal,ob,R,model,evalParam,path);

if isempty(evalDB)
    disp('no path to goal!!');
    u=[0;0];return;
end

%各評価関数の正規化
evalDB=NormalizeEval(evalDB);

%最終評価値の計算
feval=[];
for id=1:length(evalDB(:,1))
    feval=[feval;evalParam(1:3)*evalDB(id,3:5)'];
end
evalDB=[evalDB feval];

[maxv,ind]=max(feval);%最も評価値が大きい入力値のインデックスを計算
u=evalDB(ind,1:2)';%評価値が高い入力値を返す
end


function [evalDB,trajDB]=Evaluation(x,Vr,goal,ob,R,model,evalParam,path)
%各パスに対して評価値を計算する関数
evalDB=[];
trajDB=[];

for vt=Vr(1):model(5):Vr(2)
    for ot=Vr(3):model(6):Vr(4)
        %軌跡の推定
        [xt,traj]=GenerateTrajectory(x,vt,ot,evalParam(4),model);
        %各評価関数の計算
        heading=CalcHeadingEval(xt,goal);
        dist=CalcDistEval(xt,ob,R);
        path_dist=CalcPathDistEval(xt,path);
        vel=abs(vt);
        evalDB=[evalDB;[vt ot heading dist vel path_dist]];
        trajDB=[trajDB;traj];     
    end
end
end

function EvalDB=NormalizeEval(EvalDB)
%評価値を正規化する関数
if sum(EvalDB(:,3))~=0
    EvalDB(:,3)=EvalDB(:,3)/sum(EvalDB(:,3));
end
if sum(EvalDB(:,4))~=0
    EvalDB(:,4)=EvalDB(:,4)/sum(EvalDB(:,4));
end
if sum(EvalDB(:,5))~=0
    EvalDB(:,5)=EvalDB(:,5)/sum(EvalDB(:,5));
end
if sum(EvalDB(:,6))~=0
   EvalDB(:,6)=EvalDB(:,6)/sum(EvalDB(:,6));
end
end

function [x,traj]=GenerateTrajectory(x,vt,ot,evaldt,model)
%軌跡データを作成する関数
global dt;
time=0;
u=[vt;ot];%入力値
traj=x;%軌跡データ
while time<=evaldt
    time=time+dt;%シミュレーション時間の更新
    x=f(x,u);%運動モデルによる推移
    traj=[traj x];
end
end

function stopDist=CalcBreakingDist(vel,model)
%現在の速度から力学モデルに従って制動距離を計算する関数
global dt;
stopDist=0;
while vel>0
    stopDist=stopDist+vel*dt;%制動距離の計算
    vel=vel-model(3)*dt;%最高原則
end
end

function dist=CalcDistEval(x,ob,R)
%障害物との距離評価値を計算する関数

dist=2;
for io=1:length(ob(:,1))
    disttmp=norm(ob(io,:)-x(1:2)')-R;%パスの位置と障害物とのノルム誤差を計算
    if dist>disttmp%最小値を見つける
        dist=disttmp;
    end
end
end

function heading=CalcHeadingEval(x,goal)
%headingの評価関数を計算する関数

theta=toDegree(x(3));%ロボットの方位
goalTheta=toDegree(atan2(goal(2)-x(2),goal(1)-x(1)));%ゴールの方位

if goalTheta>theta
    targetTheta=goalTheta-theta;%ゴールまでの方位差分[deg]
else
    targetTheta=theta-goalTheta;%ゴールまでの方位差分[deg]
end

heading=180-targetTheta;
end

function path_dist=CalcPathDistEval(x,path)

theta=toDegree(x(3));
for i=1:length(path(:,1))
    if path(i,2)>x(2)
        Goal_path=[path(i+7,1),path(i+7,2)];
        break;
    end
end
pathTheta=toDegree(atan2(Goal_path(2)-x(2),Goal_path(1)-x(1)));

if pathTheta>theta
    targetTheta=pathTheta-theta;%ゴールまでの方位差分[deg]
else
    targetTheta=theta-pathTheta;%ゴールまでの方位差分[deg]
end

path_dist=180-targetTheta;
end

function Vr=CalcDynamicWindow(x,model)
%モデルと現在の状態からDyamicWindowを計算
global dt;
%車両モデルによるWindow
Vs=[0 model(1) -model(2) model(2)];

%運動モデルによるWindow
Vd=[x(4)-model(3)*dt x(4)+model(3)*dt x(5)-model(4)*dt x(5)+model(4)*dt];

%最終的なDynamic Windowの計算
Vtmp=[Vs;Vd];
Vr=[max(Vtmp(:,1)) min(Vtmp(:,2)) max(Vtmp(:,3)) min(Vtmp(:,4))];
%[vmin,vmax,ωmin,ωmax]
end

function x = f(x, u)
% Motion Model
global dt;
 
F = [1 0 0 0 0
     0 1 0 0 0
     0 0 1 0 0
     0 0 0 0 0
     0 0 0 0 0];
 
B = [dt*cos(x(3)) 0
    dt*sin(x(3)) 0
    0 dt
    1 0
    0 1];

x= F*x+B*u;
end

function radian = toRadian(degree)
% degree to radian
radian = degree/180*pi;
end

function degree = toDegree(radian)
% radian to degree
degree = radian/pi*180;
end

