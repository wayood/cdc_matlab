%% 軌道補正アルゴリズム
clear all;
%% 初期宣言
p.start = [0;0];
GOAL = [0;17];
global glo_obs;
global glo_gosa_obs;
global glo_obs_init;
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

[glo_gosa_obs]=gosamodel(glo_obs,p.start);%経路生成時のLM座標
gosa_plot = graph(glo_gosa_obs,p,glo_rand_size);
plot(GOAL(1,1),GOAL(2,1),'g:o','MarkerSize',10);
hold on;

%wpを手動で設定
[x,y]=ginput;
wp=[x.';
    y.'];
global wp_init;
wp = [wp GOAL];
wp_init = wp;
glo_obs_init = glo_obs;

two_init = f_twopoint(wp(:,1),p.start());
for i=1:length(wp(1,:))-1
  two(i) = f_twopoint(wp(:,i+1),wp(:,i));
end
hold on;

for j=1:length(wp(1,:))
  wp_kill(j) = plot(wp(1,j),wp(2,j),'g:o','MarkerSize',10);
  hold on;
end

%初期経路での評価
[po,sum_po] = initila_potential(glo_gosa_obs,wp,glo_rand_size);

glo_start=p.start;
%%  逐次補正プログラム
%軌道補正経路 リアルタイム

%初期化
i=1;
po_cdc=[];
slip=50;
dt=0.1;
count=1;
sr=0;
p_cdc=0;
p_cd=0;
p_init=0;
p_i=0;
v=5;
obs = ob_round(glo_gosa_obs,glo_rand_size);
ang = pi/2;

%ナビゲーション
while 1
   [wp,p.start,ang,kill_obs] = DynamicWindowApproach_for_cdc(p.start.',obs.',i,wp,ang);
   
   if i == length(wp(1,:))
       disp("Finish !!");
       break;
   end
   
   i=i+1;
   count=count+1;
   pause(0.001);
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
function [up_obs]=gosamodel(obs,start)
 for i=1:length(obs(1,:))
    [~,leng] = cart2pol(obs(1,i),obs(2,i));
    r(i)=len(obs(:,i).',start.');
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
    global glo_obs;
    global glo_obs_init;
    
    sliprate(ang_wp,v);
    up_obs(1,:)=obs(1,:)+glo_slip_x;
    up_obs(2,:)=obs(2,:)+glo_slip_y;
    glo_obs(1,:) = glo_obs_init(1,:) + glo_slip_x;
    glo_obs(2,:) = glo_obs_init(2,:) + glo_slip_y;
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

%% ポテンシャル場で初期経路を評価
function [po,sum_po] = initila_potential(glo_gosa_obs,wp,glo_rand_size)
    po = [];
    for i = 1:length(wp(1,:))-1  
        vec = wp(:,i+1) - wp(:,i);
        [ang,~] = cart2pol(vec(1,1),vec(2,1));
        x = linspace(wp(1,i),wp(1,i+1),100);
        y = linspace(wp(2,i),wp(2,i+1),100);
        for j = 1:100
            start = [x(j)
                     y(j)];
            [~,sen_num,gosa_obs] = sensor_range(glo_gosa_obs,start,ang);
            [rand_size,~] = sensor_judge(glo_gosa_obs,sen_num,glo_rand_size);
            po_st = potential(gosa_obs,start,rand_size);
            po = [po po_st];
            sum_po = sum(po)/length(po);
        end        
    end
    fprintf("Initial Potential evaluation --> %f\n",sum_po);
    currentFile = sprintf('./potential/potential.mat');
    save(currentFile,'po','sum_po');
end
%% 逐次的なアニメーション
%逐次的に軌道補正を表示
%{
function [b,w]=animation(LM_current,LM_first,wp,i,size)

    for j=i:length(wp(1,:))
        w(j)=plot(wp(1,j),wp(2,j),'g:o','MarkerSize',10);
    end
    %plot(p.start(1,1),p.start(2,1),'r:.','MarkerSize',3);
    %hold on;
    for j=1:length(LM_first(1,:))
        b(j)=en_plot_orange(LM_first(:,j).',size(j));
    end
    for j=1:length(LM_current(1,:))
        b(j)=en_plot_orange(LM_current(:,j).',size(j));
    end
end
%}
%距離を計算
function l=len(a,b)
l=norm(a-b);
end

%% 軌道補正計算
function [awp,k,mat_er,plan_er]=correction(lm_current,lm_first)
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
    [h,~]=size(lm_first);
    if h == 2
        lm_first(3,:) = 1;
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
    VTRate_spatial = corrcoef(V_current,V_first);
    CNRate_spatial = cond(Z_current)/cond(Z_first);  
    if h == 3
        [~,Z_t_1,~] = svd(LM_t_1);
        %VTRate_seque = sum(dot(V_current,V_t_1))/3;
        CNRate_seque = cond(Z_current)/cond(Z_t_1);
        Path_analysis_vir = [mat_er,plan_er,VTRate_spatial(1,2),CNRate_spatial,CNRate_seque];
        Path_analysis = [Path_analysis;Path_analysis_vir];
        file = sprintf("A_matrix.mat");
        save(file,"Path_analysis");
        fprintf("A matrix disperation --> cond %f,%f\n",CNRate_spatial,VTRate_spatial(1,2));
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

%% 二点間のプロット
function kill = f_twopoint(af_start,bef_start)
    start = [af_start bef_start];
    kill = plot(start(1,:),start(2,:),'-r','LineWidth',2);
    hold on;
end

%% アニメーション　円の作成
%円の塗りつぶし
function b=en_plot_blue(glo_obs,size)
 [x,y]=circle(glo_obs(1,1),glo_obs(1,2),size);
 b=fill(x,y,'b');
 hold on;
end

function b=en_plot_red(glo_obs,size)
 [x,y]=circle(glo_obs(1,1),glo_obs(1,2),size);
 b=fill(x,y,'r');
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


%% local plan
function [wp,start,ang,b] = DynamicWindowApproach_for_cdc(start,obstacle,wp_i,wp,ang)

    x=[start ang 0 0]';%ロボットの初期状態[x(m),y(m),yaw(Rad),v(m/s),ω(rad/s)]
    global glo_obs;
    global wp_init;
    %global N;
    %global NU;
    global glo_gosa_obs;
    global glo_rand_size;
    global drive_cdc;
    obstacleR=0.5;%衝突判定用の障害物の半径
    global po_cdc;
    Hz = 2;
    Goal_tor=0.2;
    
    %ロボットの力学モデル
    %[最高速度[m/s],最高回頭速度[rad/s],最高加減速度[m/ss],最高加減回頭速度[rad/ss],
    % 速度解像度[m/s],回頭速度解像度[rad/s]]
    Kinematic=[1.0,toRadian(20.0),0.2,toRadian(50.0),0.01,toRadian(1)];

    %評価関数のパラメータ [heading,dist,velocity,predictDT]
    evalParam=[0.1,0.2,0.1,3.0];

    %シミュレーション結果
    result.x=[];

    % Main loop
    for i=1:5000
        R = rem(i*0.5,Hz);
        if i == 1 && wp_i == 1 
            goal = wp_init(:,wp_i).';
            wp = wp_init;
            [me_gosa_obs]=gosa_hozon(glo_gosa_obs);
        else
            goal = wp(:,wp_i).';
        end

        %DWAによる入力値の計算
        if i==1
            [u,traj]=DynamicWindowApproach(x,Kinematic,goal,evalParam,obstacle,obstacleR);
            me_gosa_obs = glo_gosa_obs;
        else
            [u,traj]=DynamicWindowApproach(x,Kinematic,goal,evalParam,obs,obstacleR);
            delete(d_q);
            delete(d_g);
            delete(d_tr);            
            delete(L);
            delete(wp_plt);
        end

        x=f(x,u);%運動モデルによる移動

        %シミュレーション結果の保存
        result.x=[result.x; x'];
        start=[x(1);x(2)];
        drive_cdc=[drive_cdc start];

        % 誤差の検出と推定
        [cur_obs]=gosa_move(me_gosa_obs,start,x(3),u(1,1));
        [ang_wp,sen_num,up_obs]=sensor_range(cur_obs,start,x(3));
        [rand_size,gosa_obs]=sensor_judge(glo_gosa_obs,sen_num,glo_rand_size);
        glo_obs(3,:) = 1;
        [~,glo_range_obs]=sensor_judge(glo_obs,sen_num,glo_rand_size);
        glo_obs(3,:) = [];
        glo_range_obs(3,:) = [];
        [up_obs]=gosamodel(glo_range_obs,start);
        up_obs(3,:) = [];
        
        % 障害物の座標格納
        ob=ob_round(up_obs,rand_size);
        obs=ob.';

        %LMの変化量を図示
        L = lm_line(up_obs,gosa_obs);

        %ポテンシャル遷移で評価
        
        po_cdc_st=potential(up_obs,start,rand_size);
        po_cdc = [po_cdc po_cdc_st];
        sum_po_cdc=sum(po_cdc)/length(po_cdc);
        currentFile = sprintf('./potential/potential_cdc.mat');
        save(currentFile,'glo_obs','glo_gosa_obs','glo_rand_size','drive_cdc','po_cdc','sum_po_cdc');

        if R == 0
            [wp,k,mat_er,plan_er]=correction(up_obs,gosa_obs);
            if flag == 1
                ang = x(3);
                %プロットポイントコメントアウト部分
                delete(wp_plt);
                delete(L);
                delete(r);
                break;
            end
        end
        for j=wp_i:length(wp(1,:))
          wp_plt(j) = plot(wp(1,j),wp(2,j),'g:o','MarkerSize',10);
          hold on;
        end

        %プロットポイントコメントアウト部分


        if i>1
            delete(b);
            delete(r);
        end
        for obs_i = 1:length(glo_obs(1,:))
            r(obs_i)=en_plot_red(glo_obs(:,obs_i).',glo_rand_size(obs_i));
        end
        for j=1:length(up_obs(1,:))
            b(j)=en_plot_orange(up_obs(:,j).',rand_size(j));
        end

        %ゴール判定
        if norm(x(1:2)-goal')<Goal_tor
            disp('Arrive Goal!!');
            ang = x(3);
            %プロットポイントコメントアウト部分
            delete(wp_plt);
            delete(L);
            delete(r);
            if length(wp(1,:)) == wp_i
                break;
            end
            delete(b);        
            break;
        end


        if  i > 30 && abs(result.x(length(result.x(:,1)),1) - result.x(length(result.x(:,1))-30,1)) < 1.0 && abs(result.x(length(result.x(:,1)),2) - result.x(length(result.x(:,1))-30,2)) < 1.0 
            obstacleR=0.1;
        end

        if  i > 40 && abs(result.x(length(result.x(:,1)),1) - result.x(length(result.x(:,1))-40,1)) < 1.0 && abs(result.x(length(result.x(:,1)),2) - result.x(length(result.x(:,1))-40,2)) < 1.0 
            obstacleR=0.0;
            evalParam(6)=0.3;
        end

        if  i > 100 && abs(result.x(length(result.x(:,1)),1) - result.x(length(result.x(:,1))-100,1)) < 1.0 && abs(result.x(length(result.x(:,1)),2) - result.x(length(result.x(:,1))-100,2)) < 1.0 
            Goal_tor=5.0;
        end

        if  i > 150 && abs(result.x(length(result.x(:,1)),1) - result.x(length(result.x(:,1))-150,1)) < 1.0 && abs(result.x(length(result.x(:,1)),2) - result.x(length(result.x(:,1))-150,2)) < 1.0 
            disp('Skip Waypoint');
            %プロットポイントコメントアウト部分
               delete(b);
            break;
        end

        %プロットポイントコメントアウト部分

        if i>1    
            delete(d_x);   
        end

            %====Animation====
         ArrowLength=0.5;%矢印の長さ
         %ロボット
         d_q=quiver(x(1),x(2),ArrowLength*cos(x(3)),ArrowLength*sin(x(3)),'ok');
         hold on;
         d_x=plot(result.x(:,1),result.x(:,2),'-b','LineWidth',2);
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
         %drawnow;
         pause(0.01);
    end

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