%% 軌道補正アルゴリズム
clear all;
numFiles=30;
global N;
for N=1:numFiles
hold off;

global cdc_length;
cdc_length=0.5;%m単位で補正ポイントを設定
global glo_obs;
global glo_gosa_obs;
global glo_rand_size;
%% パラメータ設定
p.start=[0;0];
p.goal=[0;50];
for i=1:700
    ran_x=-50+100*rand;
    ran_y=300*rand;
    glo_rand_size(i)=0.3+1.4*rand;
    glo_obs(1,i)=ran_x;
    glo_obs(2,i)=ran_y;
end

[glo_gosa_obs]=gosamodel(glo_obs,p);%経路生成時のLM座標
obs=ob_round(glo_gosa_obs,glo_rand_size);
%graph(glo_gosa_obs,p,glo_rand_size);
DynamicWindowApproachSample_k(p.start.',p.goal.',obs.')
%{
hold on;
for io=1:length(glo_obs(1,:))
  en_plot_red(glo_obs(:,io).',glo_rand_size(io));
  hold on;
end
%}
disp("finish !!");
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


%距離を計算
function l=len(a,b)
l=norm(a-b);
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
    ylim([0 200]);
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

%% 経路生成時のLM座標の視野角判定
function [rand_size,cur_gosa_obs]=sensor_judge(gosa_obs,sen_num,glo_rand_size)
   for i=1:length(gosa_obs(1,:))
      if sen_num(i)==1
          gosa_obs(:,i)=[-1;-1;-1];
          glo_rand_size(i)=0;
      end
   end
    idx = gosa_obs(1,:)==-1 & gosa_obs(2,:) == -1 & gosa_obs(3,:) == -1;
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

%円の関数
function [r_x,r_y]=circle(x,y,r)
 t=linspace(0,2*pi,100);
 r_x=r*cos(t)+x;
 r_y=r*sin(t)+y;
end

%% 障害物の円座標を格納
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
function DynamicWindowApproachSample_k(start,goal,obstacle)

x=[start pi/2 0 0]';%ロボットの初期状態[x(m),y(m),yaw(Rad),v(m/s),ω(rad/s)]
global glo_obs;
global glo_gosa_obs;
global glo_rand_size;
global drive_cdc;
obstacleR=0.2;%衝突判定用の障害物の半径
global dt; 
global N;
dt=0.1;%刻み時間[s]
drive_cdc=[];
i_po=1;

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
[u,traj]=DynamicWindowApproach(x,Kinematic,goal,evalParam,obstacle,obstacleR);
x=f(x,u);%運動モデルによる移動

%シミュレーション結果の保存
result.x=[result.x; x'];
start=[x(1),x(2)];
s_x=[x(1);x(2)];
drive_cdc=[drive_cdc s_x];
[~,sen_num,gosa_obs]=sensor_range(glo_gosa_obs,start.',goal.');
[rand_size,~]=sensor_judge(glo_gosa_obs,sen_num,glo_rand_size);
po_cruise(i_po)=potential(gosa_obs,s_x,rand_size);
sum_po_cruise=sum(po_cruise)/i_po;
i_po=i_po+1;

%% saveする際のファイル名設定
currentFile = sprintf('potential_cruise_%d_50.mat',N);
save(currentFile,'po_cruise','sum_po_cruise');
currentFile = sprintf('path_interp_%d_50.mat',N);
save(currentFile,'drive_cdc','glo_obs','glo_gosa_obs','glo_rand_size');

if i>20
    
    V_x=var(drive_cdc(1,[i-19 i]));
    V_y=var(drive_cdc(2,[i-19 i]));

    if V_x < 0.1 && V_y <0.1
        disp("Fuck");
        break;
    end
end
%{
if i>1
    delete(d_q);
    delete(d_g);
    delete(d_tr);
end
%}


%ゴール判定
if norm(x(1:2)-goal')<1.0
    disp('Arrive Goal!!');
    s=[x(1);x(2)];
    break;
end

%{
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
%}

 drawnow;
end
%toc
end
%movie2avi(mov,'movie.avi');
 

function [u,trajDB]=DynamicWindowApproach(x,model,goal,evalParam,ob,R)
%DWAによる入力値の計算をする関数

%Dynamic Window[vmin,vmax,ωmin,ωmax]の作成
Vr=CalcDynamicWindow(x,model);
%評価関数の計算
[evalDB,trajDB]=Evaluation(x,Vr,goal,ob,R,model,evalParam);

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


function [evalDB,trajDB]=Evaluation(x,Vr,goal,ob,R,model,evalParam)
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
        vel=abs(vt);
        evalDB=[evalDB;[vt ot heading dist vel]];
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
for i=1:length(path(1,:))
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

