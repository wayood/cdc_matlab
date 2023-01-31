%% 軌道補正アルゴリズム
clear all;
%try
fig = figure;
fig.Color = 'white';


%% 初期宣言
p.start = [0;0];
global range_base;
global glo_obs;
global glo_gosa_obs;
global glo_obs_init;
global glo_rand_size;
global glo_add_obs_init;
global add_obs_init;
global ground_add_obs;
global add_obs_rand_size;
global drive_cdc;
global dt;
global slip;
global glo_slip_x;
global glo_slip_y;
global po_cdc;
global lm_cur_1;
global Path_analysis;
global flag_add;
global ax;
global ButtonState;
global add_path;
global obs_list;
global flag_wp_continue;
global wp_add;
global lm_add_stock;
global lm_first_stock;
global add_count;
global lm_add_current_stock;
global add_obs_size;
global first_obs_size;
global lm_first_current_stock;
global add_slip;
add_slip = {};
add_count = 0;
lm_first_stock = {};
lm_add_stock = {};
wp_add = {};
lm_add_current_stock = {};
add_obs_size = {};
first_obs_size = {};
lm_first_current_stock = {};
flag_wp_continue = 0;
add_path.State = 'wait';
flag_add = 0;
Path_analysis =[];
lm_cur_1 = [];
glo_slip_x = 0;
glo_slip_y = 0;
drive_cdc=[];
ang_wp = pi/2;
range_base=10;
i=1;
j=1;
ax = axes(fig);
ax.GridLineStyle = '-';
ax.MinorGridAlpha = 0.1;
ax.LineWidth = 1.5;
ax.FontSize= 14;
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';

fig.WindowButtonDownFcn = @WindowButtonDownFcn_fig;
fig.WindowButtonUpFcn = @WindowButtonUpFcn_fig;


%% 障害物
while 1
    ran_x=-20+40*rand;
    ran_y=100*rand;
    [ang,l] = cart2pol(ran_x,ran_y);   
    
    if l <= range_base  && (7*pi)/36 <= ang && ang <= (29*pi)/36
        glo_rand_size(i)=0.005+0.5*rand;
        glo_obs(1,i)=ran_x;
        glo_obs(2,i)=ran_y;
        i=i+1;
    elseif (7*pi)/36 <= ang && ang <= (29*pi)/36 && ran_y < range_base + 20
        add_obs_rand_size(j) = 0.005+0.5*rand;
        ground_add_obs(:,j) = [ran_x;ran_y];
        j=j+1;
    end

    if i >= 30 && j >= 75
        break;
    end
end

[glo_gosa_obs]=gosamodel(glo_obs,p.start,ang_wp);%経路生成時のLM座標
[add_obs_init]=gosamodel(ground_add_obs,p.start,ang_wp);
glo_add_obs_init = ground_add_obs;
glo_obs_init = glo_obs;

gosa_plot = graph_first(glo_gosa_obs,p,glo_rand_size);


%wpを手動で設定
[x,y]=ginput;
wp=[x.';
    y.'];
global wp_init;
wp_init = wp;

currentFile = sprintf('./wp/path_obstacle.mat');
save(currentFile,'glo_obs_init','wp_init','glo_gosa_obs','glo_rand_size');

two_init = f_twopoint(wp(:,1),p.start());
for i=1:length(wp(1,:))-1
  two(i) = f_twopoint(wp(:,i+1),wp(:,i));
end
hold on;

%{
for j=1:length(wp(1,:))
  wp_kill(j) = plot(wp(1,j),wp(2,j),'g:o','MarkerSize',10);
  hold on;
end
%}

%初期経路での評価
[po,sum_po] = initial_potential(glo_gosa_obs,wp,glo_rand_size);

glo_start=p.start;
%%  逐次補正プログラム
%軌道補正経路 リアルタイム

%初期化
i=1;
po_cdc=[];
slip = 30;
dt = 0.1;
count = 1;
elapsed_count = 1;

obs = ob_round(glo_gosa_obs,glo_rand_size);
ang = pi/2;
wp_add_count = 1;
wp_add_array(1).wp = [];
wp_add_array(1).count = [];
wp_add_array(1).lm_add = [];
wp_add_array(1).lm_add_range = [];
ButtonState = false;
obs_list = obs.';

%ナビゲーション
while 1
   [wp,p.start,ang,flag,kill_obs,current_obs,current_rand_size] = DynamicWindowApproach_for_cdc(p.start.',obs_list,i,wp,ang,wp_add_array); 
   if i == length(wp(1,:))
        disp("Finish");
        break;
   end
   
   %{
   [x,y]=ginput;
    wp_s=[x.';
        y.'];
    wp = [wp wp_s];
   %}

   %add_path = parfeval(backgroundPool,@DynamicWindowApproach_global,2,wp(:,end-1).',wp(:,end).',obs);
    
    
   if flag == 1
       disp("Please add waypoint !!");
       
       delete(gosa_plot);
       delete(two);
       delete(two_init);
       %delete(wp_kill);
       
       
       tic
       [rrt_path] =  RRTstar3D(p.start,wp(:,i),current_obs,current_rand_size);
       elapsedTime_astar(elapsed_count) = toc;
       plot([p.start(1,1) ;rrt_path(:,1) ;wp(1,i)],[p.start(2,1) ;rrt_path(:,2);wp(2,i)],"-m",'LineWidth',2);
       hold on;
       
       [po_rrt_star,sum_po_rrt_star] = initial_potential(current_obs,[p.start rrt_path(:,1:2).' wp(:,i)],current_rand_size);
       R_rrt_star = correlation([p.start wp(:,i)],[p.start rrt_path(:,1:2).' wp(:,i)]);
       fprintf("RRT star : time %f [sec] , safe %f ,correlation %f [%]\n\n",elapsedTime_astar(1),sum_po_rrt_star,R_rrt_star(1,2));
       
       tic
       [wp_replan] = voronoi_waypoint_generation(p.start,wp(:,i),current_obs(1:2,:));
       elapsedTime_add(elapsed_count) = toc;
       elapsed_count = elapsed_count + 1;
       [po_voronoi,sum_po_voronoi] = initial_potential(current_obs,[p.start wp_replan],current_rand_size);
       R_voronoi = correlation([p.start wp(:,i)],[p.start wp_replan]);
       fprintf("Voronoi graph : time %f [sec] , safe %f ,correlation %f [%]\n\n",elapsedTime_add(1),sum_po_voronoi,R_voronoi(1,2));
       
       % [po_st,sum_po_st] = initila_potential(glo_gosa_obs,wp,glo_rand_size);
       pl_wp=[wp_init(:,1:i-1) wp_replan wp_init(:,i+2:end)];
       two = plot([wp_replan(1,:) wp_init(1,i+2:end)],[wp_replan(2,:) wp_init(2,i+2:end)],'-r','LineWidth',2);
       wp_kill = plot(pl_wp(1,:),pl_wp(2,:),'g:o','MarkerSize',10);
       wp_replan = wp_replan(:,2:end);
       wp = [wp_init(:,1:i-1) wp_replan wp_init(:,i+2:end)];
       
       wp_add_array(wp_add_count).wp = wp_replan;
       wp_add_array(wp_add_count).count = i;
       lm_add_array(1,:) = glo_gosa_obs(1,:);
       lm_add_array(2,:) = glo_gosa_obs(2,:);
       wp_add_array(wp_add_count).lm_add = lm_add_array;
       
       i = i-1;
       
       obs = ob_round(glo_gosa_obs,glo_rand_size);
       delete(kill_obs);
     
       for plt_cou = 1:length(glo_gosa_obs(1,:))
          [x,y]=circle(glo_gosa_obs(1,plt_cou),glo_gosa_obs(2,plt_cou),glo_rand_size(plt_cou));
          gosa_plot(plt_cou)=fill(x,y,'b','FaceAlpha',.3,'EdgeAlpha',.3);
          hold on;
       end
       
       wp_add_count = wp_add_count + 1;
       wp_init = wp;
       
   end
   i=i+1;
   count=count+1;

end
%{
catch
    cancel(add_path);
end
%}
%% callback関数の設定
function WindowButtonDownFcn_fig(~,~)

    global ax;
    global ButtonState;
    global obs_list;
    global add_path;
    global wp_init;
    global flag_wp_continue;
    global add_count;
    global flag_add;
    global wp_add;

    ButtonState = true; % 押下状態を保存
    if isempty(wp_add) == 0
        add_path = parfeval(backgroundPool,@DynamicWindowApproach_global,2,wp_add{end,:}(:,end).',[ax.CurrentPoint(1,1),ax.CurrentPoint(1,2)],obs_list);
    else
        add_path = parfeval(backgroundPool,@DynamicWindowApproach_global,2,wp_init(:,end).',[ax.CurrentPoint(1,1),ax.CurrentPoint(1,2)],obs_list);
    end
    flag_add = 1;
    if flag_wp_continue == 1
        flag_wp_continue = 0;
    end
    add_count = add_count + 1;
    plot(ax.CurrentPoint(1,1),ax.CurrentPoint(1,2),'g:o','MarkerSize',10);
    hold on;

end

function WindowButtonUpFcn_fig(~,~)
    global ButtonState;
    ButtonState = false;

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
    for i = 1:length(LM_first(1,:))
        LM = [LM_current(:,i) LM_first(:,i)];
        kill(i)=plot(LM(1,:),LM(2,:),'-g','LineWidth',1.5);
        hold on;
    end
end
%% 誤差モデル計算　
%距離による関係性を考えた。
%これはDEMのステレオデータによる誤差を主に考えた。
function [up_obs]=gosamodel(obs,start,ang)
 global range_base;
 for i=1:length(obs(1,:))
    [~,leng] = cart2pol(obs(1,i),obs(2,i));
    r(i)=len(obs(:,i).',start.');
    if r(i)<1
        r(i)=0;
    end
    x_ran=-0.01+0.02*rand; %キャリブレーションや分解能での誤差を考える
    y_ran=-0.01+0.02*rand;
    if leng >= range_base
        l=18.5*10^-2*range_base^2*0.01;
        [x_error,y_error] = pol2cart(ang,l);
        up_obs(1,i)=obs(1,i)+x_error+x_ran;
        up_obs(2,i)=obs(2,i)+y_error+y_ran;
        continue;
    end
    l=18.5*10^-2*r(i)^2*0.01;
    [x_error,y_error] = pol2cart(ang,l);
    up_obs(1,i)=obs(1,i)+x_error+x_ran;
    up_obs(2,i)=obs(2,i)+y_error+y_ran;
 end
end

%誤差をロボット進行方向に考えて表示
function [up_obs]=gosa_move(obs,start,ang_wp,v)
    global glo_slip_x;
    global glo_slip_y;
    
    sliprate(ang_wp,v);
    up_obs(1,:)=obs(1,:)+glo_slip_x;
    up_obs(2,:)=obs(2,:)+glo_slip_y;

end

function R = correlation(one_wp,another_wp)
    N = length(another_wp(1,:)) - 1;
    one_x = linspace(one_wp(1,1),one_wp(1,end),N*100);
    one_y = linspace(one_wp(2,1),one_wp(2,end),N*100);
    one = [one_x
           one_y];
       
    another_x = [];
    another_y = [];
    for i = 2:length(another_wp(1,:))
        stock_x = linspace(another_wp(1,i-1),another_wp(1,i),100);
        stock_y = linspace(another_wp(2,i-1),another_wp(2,i),100);
        another_x = [another_x stock_x];
        another_y = [another_y stock_y];
    end
    another = [another_x
               another_y];
    R = corrcoef(one,another);
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
 glo_slip_x =  glo_slip_x - x;
 glo_slip_y =  glo_slip_y - y;
end

%% ポテンシャル場で初期経路を評価
function [po,sum_po] = initial_potential(glo_gosa_obs,wp,glo_rand_size)
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
function [wp_new,k,mat_er,plan_er,flag]=compensation(lm_current,lm_first,lm_add_current,lm_add_init)
    global wp_init;
    global lm_cur_1;
    global wp_add;
    global A_n;
    wp_init(3,:)=1;

    if isempty(lm_current) == 0
        lm_current(3,:) = 1;
        lm_first(3,:) = 1; 
        A=lm_current*pinv(lm_first);
        flag = 0;
        if isempty(A_n)==0
         [mat_er,plan_er,flag] = A_matrix(A,lm_current,lm_first,A_n,lm_cur_1);
        else
            mat_er=0;
            plan_er=0;
        end
        lm_cur_1 = lm_current;
        A_n=A;
        wp_new=A*wp_init;
        wp_init(3,:) = [];
        k=cond(A,2);
    else
        wp_new = wp_init;
        wp_new(3,:)=[];
        k = 0;
        mat_er = 0;
        plan_er = 0;
        flag = 0;
    end
    
    if isempty(wp_add) == 0
        for i = 1:length(wp_add)
            wp_add{i,:}(3,:) = 1;
            
            lm_add_current{i}(3,:) = 1;
            lm_add_init{i}(3,:) = 1;
            
            A_add = lm_add_current{i}*pinv(lm_add_init{i});
            wp_add_new = A_add*wp_add{i,:};
            wp_new = [wp_new wp_add_new];
            wp_add{i,:}(3,:) = [];
        end
        wp_new(3,:)=[];
    end
    
end

%% A行列解析
%カリタニさん参照

function [matrix_error,matrix_timespace_error,flag]=A_matrix(A,LM_current,LM_first,A_n,LM_t_1)
    global Path_analysis;
    matrix_error = cond(A)*norm(A*LM_first-LM_current)/norm(A*LM_first);
    matrix_timespace_error = cond(A_n)*norm(A_n-A)/norm(A_n);
    [h,~] = size(LM_t_1);
    [~,S_first,V_first] = svd(LM_first);
    [~,S_current,V_current] = svd(LM_current);
    VTRate_spatial = corrcoef(V_current,V_first);
    CNRate_spatial = cond(S_current)/cond(S_first);
    [U_distortion,S_distortion,V_distortion] = svd(A);
    
    rotation = U_distortion*V_distortion;
    
    theta_z = atan(-rotation(1,2)/rotation(1,1));
    parallel_y = 2*(-rotation(3,2)/rotation(3,3));
    parallel_x = 2*((rotation(3,1)+rotation(3,2))/(rotation(2,1)+rotation(2,2)));
    
    rotation_theta_1 = asin(U_distortion(1,1));
    rotation_theta_2 = asin(V_distortion(1,1));
    streching_x = S_distortion(1,1);
    streching_y = S_distortion(2,2);
    flag = 0;
    
    if matrix_error > 1.0 || matrix_timespace_error > 0.3
        flag =0;
    end
    
    if h == 3
        [~,Z_t_1,~] = svd(LM_t_1);
        %VTRate_seque = sum(dot(V_current,V_t_1))/3;
        CNRate_seque = cond(S_current)/cond(Z_t_1);
        Path_analysis_vir = [matrix_error,matrix_timespace_error,parallel_x,parallel_y,theta_z,VTRate_spatial(1,2),CNRate_spatial,CNRate_seque,cond(A)];
        Path_analysis = [Path_analysis;Path_analysis_vir];
        file = sprintf("A_matrix.mat");
        save(file,"Path_analysis");
        fprintf("A matrix disperation (Utsuno proposal)\n  --> Error Size %f, Error direction %f\n\n",CNRate_spatial,VTRate_spatial(1,2));
        fprintf("A matrix disperation (Karitani proposal)\n  --> Matrix Error %f, Timespace Error %f, Parallel Translation(x) %f, Parallel Translation(y) %f,Rotation %f \n\n",matrix_error,matrix_timespace_error,parallel_x,parallel_y,theta_z);
        fprintf("A matrix disperation (main proposed)\n  --> Rotation %f, Streching x %f, Streching Y %f\n\n",rotation_theta_1+rotation_theta_2,streching_x,streching_y);
    end
end

%% センサ検知前位置のグラフを表示
function [b] = graph_first(glo_gosa_obs,p,size)
    
    for i=1:length(glo_gosa_obs(1,:))
        b(i) = en_plot_blue(glo_gosa_obs(:,i).',size(i));
    end
    plot(p.start(1,1),p.start(2,1),'b:.','MarkerSize',5);
    hold on;
    grid on;
    xlabel('$x[m]$', 'Interpreter', 'latex');
    ylabel('$y[m]$', 'Interpreter', 'latex');
    xlim([-20 20]);
    ylim([0 30]);
end


%% 経路生成時のLM座標の視野角判定
function [rand_size,cur_gosa_obs]=sensor_judge(gosa_obs,sen_num,glo_rand_size)
   for i=1:length(gosa_obs(1,:))
      if sen_num(i)==1
          gosa_obs(:,i)=[-1;-1];
          glo_rand_size(i)=0;
      end
   end
    idx = gosa_obs(1,:)==-1 & gosa_obs(2,:) == -1;
    idy = glo_rand_size(1,:)==0;
    rand_size=glo_rand_size(~idy);
    cur_gosa_obs=gosa_obs(:,~idx);
end

%% 視野角を考慮
function [ang,sen_num,sen_obs]=sensor_range(obs,start,ang)
   %視野の射程距離
   global range_base;
    [rows,~] = size(obs);
   if rows == 3
       obs(3,:)=[];
   end
   if ang > pi
       ang_pol = ang - 2*pi;
   elseif ang > 2*pi
       r = rem(ang,2*pi);
       if r > pi
           ang_pol = r - 2*pi;
       else
           ang_pol = r;
       end
   else
       ang_pol = ang;
   end
   range_min = ang_pol-(11*pi/36);
   range_max = ang_pol+(11*pi/36);
   if abs(range_min) > pi && ang_pol < 0
       range_plus_min = 2*pi - abs(range_min);
       range_plus_max = pi;
       range_minus_max = range_max;
       range_minus_min = -pi;    
   elseif range_max > pi && ang_pol > 0
       range_plus_min = range_min;
       range_plus_max = pi;
       range_minus_max = range_max - 2*pi;
       range_minus_min = -pi;
   else
       range_plus_min = range_min;
       range_plus_max = range_max;
       range_minus_max = -2;
       range_minus_min = -1;
   end
   
   for i=1:length(obs(1,:))
        [ang_obs,range_length] = cart2pol(obs(1,i)-start(1,1),obs(2,i)-start(2,1)); 
        if ang_obs>range_plus_min && ang_obs<range_plus_max && range_length<range_base
            sen_num(i)=0;
        elseif ang_obs>range_minus_min && ang_obs<range_minus_max && range_length<range_base
            sen_num(i)=0;
        else
            obs(:,i)=[-1;-1];
            sen_num(i)=1;
        end
   end
      idx = obs(1,:)== -1 & obs(2,:) == -1;
      sen_obs = obs(:,~idx);
end

function kill = plot_sensor_range(start,ang)
    
    global range_base;
    range_minus_min = 0;
   if ang > pi
       ang_pol = ang - 2*pi;
   elseif ang > 2*pi
       r = rem(ang,2*pi);
       if r > pi
           ang_pol = r - 2*pi;
       else
           ang_pol = r;
       end
   else
       ang_pol = ang;
   end
   range_min = ang_pol-(11*pi/36);
   range_max = ang_pol+(11*pi/36);
   if abs(range_min) > pi && ang_pol < 0
       range_plus_min = 2*pi - abs(range_min);
       range_plus_max = pi;
       range_minus_max = range_max;
       range_minus_min = -pi;    
   elseif range_max > pi && ang_pol > 0
       range_plus_min = range_min;
       range_plus_max = pi;
       range_minus_max = range_max - 2*pi;
       range_minus_min = -pi;
   else
       range_plus_min = range_min;
       range_plus_max = range_max;
       side_range = [start(1,1),start(2,1)
                     range_base*cos(range_min) + start(1,1),range_base*sin(range_min) + start(2,1)
                     range_base*cos(range_max) + start(1,1),range_base*sin(range_max) + start(2,1)];
   end

   if isempty(find(range_minus_min == 0)) == 0
       theta = linspace(range_min,range_max);
       kill(1) = plot(range_base*cos(theta) + start(1,1),range_base*sin(theta)+ start(2,1),'Color',[0,0,0],'LineWidth',1);
       hold on;
   else
       theta_plus = linspace(range_plus_min,range_plus_max);
       theta_minus = linspace(range_minus_min,range_minus_max);
       side_range = [start(1,1),start(2,1)
                     range_base*cos(range_plus_min) + start(1,1),range_base*sin(range_plus_min) + start(2,1)
                     range_base*cos(range_minus_max) + start(1,1),range_base*sin(range_minus_max) + start(2,1)];
       kill(1) = plot(range_base*cos(theta_minus) + start(1,1),range_base*sin(theta_minus)+ start(2,1),'Color',[0,0,0],'LineWidth',1);
       hold on;
       kill(2) = plot(range_base*cos(theta_plus) + start(1,1),range_base*sin(theta_plus)+ start(2,1),'Color',[0,0,0],'LineWidth',1);
       hold on;
   end
   kill(end+1) = plot(side_range(1:2,1),side_range(1:2,2),'Color',[0,0,0],'LineWidth',1);
    hold on;
    side_range(2,:)=[];
    kill(end+1) = plot(side_range(1:2,1),side_range(1:2,2),'Color',[0,0,0],'LineWidth',1);
    hold on;
end
%% 二点間のプロット
function kill = f_twopoint(af_start,bef_start)
    start = [af_start bef_start];
    kill = plot(start(1,:),start(2,:),'Color',[1, 0, 0, 0.2],'LineWidth',3);
    hold on;
end

%% アニメーション　円の作成
%円の塗りつぶし
function b=circle_plot(obstacle,size,color)
 [x,y]=circle(obstacle(1,1),obstacle(1,2),size);
 b=fill(x,y,color);
 hold on;
end

function b=en_plot_blue(glo_obs,size)
 [x,y]=circle(glo_obs(1,1),glo_obs(1,2),size);
 b=fill(x,y,'b','FaceAlpha',.3,'EdgeAlpha',.3);
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
 
 for i=1:length(obs(1,:))
    up_obs(1,i)=obs(1,i)+glo_slip_x;
    up_obs(2,i)=obs(2,i)+glo_slip_y;
 end
 
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
function [wp,start,ang,flag,b,up_obs,rand_size] = DynamicWindowApproach_for_cdc(start,obstacle,wp_i,wp,ang,wp_add_array)

    x=[start ang 0 0]';%ロボットの初期状態[x(m),y(m),yaw(Rad),v(m/s),ω(rad/s)]
    global glo_obs;
    global wp_init;
    global add_path;
    global glo_gosa_obs;
    global glo_rand_size;
    global drive_cdc;
    global po_cdc;
    global obs_list;
    global flag_wp_continue;
    global add_count;
    global glo_add_obs_init;
    global ground_add_obs;
    global add_obs_rand_size;
    global add_obs_init;
    global flag_add;
    global wp_add;
    global lm_add_stock;
    global lm_first_stock;
    global glo_slip_x;
    global glo_slip_y;
    global glo_obs_init;
    global lm_add_current_stock;
    global add_obs_size;
    global first_obs_size;
    global lm_first_current_stock;
    global add_slip;
    Hz = 2;
    Goal_tor = 0.5;
    flag = 0;
    obstacleR=0.5;%衝突判定用の障害物の半径
    mintblue = [0.3010 0.7450 0.9330];
    add_obs = [];
    add_flag = 0;
    add_obs_plot = [];

    %ロボットの力学モデル
    %[最高速度[m/s],最高回頭速度[rad/s],最高加減速度[m/ss],最高加減回頭速度[rad/ss],
    % 速度解像度[m/s],回頭速度解像度[rad/s]]
    Kinematic=[1.0,toRadian(20.0),0.2,toRadian(50.0),0.01,toRadian(1)];

    %評価関数のパラメータ [heading,dist,velocity,predictDT]
    evalParam=[0.5,0.0,0.5,2.0];

    %シミュレーション結果
    result.x=[];
    select_add_num = 1;

    % Main loop
    for i=1:5000
        first_gosa_obs = {};
        current_init_obs = {};
        add_gosa_obs = {};
        current_add_obs = {};
        lm_add_init = {};
        lm_add_current = {};
        lm_first_size = {};
        lm_current_size ={};
        lm_add_size = {};

        R = rem(i*0.5,Hz);
        if i == 1 && wp_i == 1 
            goal = wp_init(:,wp_i).';
            wp = wp_init;
            [obs_init_stock]=gosa_hozon(glo_gosa_obs);
            [add_obs_init_stock]=gosa_hozon(add_obs_init);
        else
            goal = wp(1:2,wp_i).';
        end

        
        %DWAによる入力値の計算
        if i==1
            obs_list = obstacle;
            [u,traj]=DynamicWindowApproach(x,Kinematic,goal,evalParam,obs_list,obstacleR);
            obs_init_stock = glo_gosa_obs;
            add_obs_init_stock = add_obs_init;
        else
            [u,traj]=DynamicWindowApproach(x,Kinematic,goal,evalParam,obs_list,obstacleR);
            delete(d_q);
            delete(d_g);
            delete(d_tr);            
            delete(L);
            delete(b);
            delete(r);
            delete(wp_plt);
            delete(kill_sensor_range);
        end
        
        if i>1 & isempty(find(select_add_num == 0)) == 0
            delete(d_obs);
        end

        x=f(x,u);%運動モデルによる移動

        %シミュレーション結果の保存
        result.x=[result.x; x'];
        start=[x(1);x(2)];
        drive_cdc=[drive_cdc start];

        % 誤差の検出と推定
        [cur_obs]=gosa_move(obs_init_stock,start,x(3),u(1,1));
        [cur_add_obs]=gosa_move(add_obs_init_stock,start,x(3),u(1,1));
        glo_obs(1,:) = glo_obs_init(1,:) + glo_slip_x;
        glo_obs(2,:) = glo_obs_init(2,:) + glo_slip_y;
        ground_add_obs(1,:) = glo_add_obs_init(1,:) + glo_slip_x;
        ground_add_obs(2,:) = glo_add_obs_init(2,:) + glo_slip_y;
        
        [ang_wp,sen_num,~]=sensor_range(cur_obs,start,x(3));
        [~,select_add_num,~]=sensor_range(cur_add_obs,start,x(3));

        [rand_size,gosa_obs]=sensor_judge(glo_gosa_obs,sen_num,glo_rand_size);

            %add時のlm初期状態と追加状態で登録
        if flag_add == 1 && strcmp(add_path.State,'running')  
            if find(sen_num) ~= 1
                [add_size,glo_range_obs]=sensor_judge(glo_obs,sen_num,glo_rand_size);
                [cur_obs]=gosamodel(glo_range_obs,start,ang_wp);
                lm_first_stock(add_count,:) = {cur_obs};
                lm_first_current_stock(add_count,:) = {cur_obs};
                first_obs_size(add_count,:) = {add_size};
            end
        end

        if flag_add == 1 & isempty(find(select_add_num == 0)) == 0  
            if strcmp(add_path.State,'running')
                [add_size,add_range_ground_obs]=sensor_judge(ground_add_obs,select_add_num,add_obs_rand_size);
                [cur_add_obs]=gosamodel(add_range_ground_obs,start,ang_wp);
                lm_add_stock(add_count,:) = {cur_add_obs};
                lm_add_current_stock(add_count,:) = {cur_add_obs};
                add_obs_size(add_count,:) = {add_size};
            end
        end
        
        if flag_add == 1
            flag_add = 0;
            add_slip(add_count) = {[glo_slip_x,glo_slip_y]}; 
            if length(lm_first_stock) ~= add_count
                lm_first_stock(add_count,:) = {zeros(2,length(glo_gosa_obs(1,:)))};
                lm_first_current_stock(add_count,:) = {zeros(2,length(glo_gosa_obs(1,:)))};
            end
            if length(lm_add_stock) ~= add_count
                lm_add_stock(add_count,:) = {zeros(2,length(add_obs_init(1,:)))};
                lm_add_current_stock(add_count,:) = {zeros(2,length(add_obs_init(1,:)))};
            end
        end

        if isempty(lm_add_stock) == 0
            for add_i = 1:add_count
                lm_add_current_stock{add_i}(1,:) = lm_add_stock{add_i}(1,:) + glo_slip_x - add_slip{add_i}(1,1);
                lm_add_current_stock{add_i}(2,:) = lm_add_stock{add_i}(2,:) + glo_slip_y - add_slip{add_i}(1,2);
                if isempty(find(cell2mat(lm_add_stock(add_i,:)),1)) == 1
                    gosa_add_obs = zeros(2,length(add_obs_init(1,:)));
                    add_obs = zeros(2,length(add_obs_init(1,:)));
                    add_range_size = zeros(length(add_obs_init(1,:)));
                else
                    [~,num_add,~]=sensor_range(lm_add_current_stock{add_i},start,x(3));
                    [add_range_size,gosa_add_obs]=sensor_judge(cell2mat(lm_add_stock(add_i,:)),num_add,add_obs_size{add_count});
                    [~,add_range_ground_obs]=sensor_judge(lm_add_current_stock{add_i},num_add,add_obs_size{add_count});
                    [add_obs]=gosamodel(add_range_ground_obs,start,ang_wp);
                end
                current_add_obs(add_i,:) = {add_obs};
                add_gosa_obs(add_i,:) = {gosa_add_obs};
                lm_current_size(add_i,:) = {add_range_size};
            end
        end

        if isempty(lm_first_stock) == 0
            for add_i = 1:add_count
                lm_first_current_stock{add_i}(1,:) = lm_first_stock{add_i}(1,:) + glo_slip_x - add_slip{add_i}(1,1);
                lm_first_current_stock{add_i}(2,:) = lm_first_stock{add_i}(2,:) + glo_slip_y - add_slip{add_i}(1,2);
                if isempty(find(cell2mat(lm_first_stock(add_i,:)),1)) == 1
                    gosa_add_obs = zeros(2,length(glo_gosa_obs(1,:)));
                    init_obs = zeros(2,length(glo_gosa_obs(1,:)));
                    first_range_size = zeros(length(glo_gosa_obs(1,:)));
                else
                    [~,num_first,~]=sensor_range(lm_first_current_stock{add_i},start,x(3));
                    [first_range_size,gosa_add_obs]=sensor_judge(cell2mat(lm_first_stock(add_i,:)),num_first,first_obs_size{add_i});
                    [~,first_range_gosa_obs]=sensor_judge(lm_first_current_stock{add_i},num_first,first_obs_size{add_i});
                    [init_obs]=gosamodel(first_range_gosa_obs,start,ang_wp);
                end
                current_init_obs(add_i,:) = {init_obs};
                first_gosa_obs(add_i,:) = {gosa_add_obs};
                lm_first_size(add_i,:) = {first_range_size};
            end
        end
        
        if isempty(lm_add_stock) == 0 && isempty(lm_first_stock) == 0   
            for add_i = 1:add_count
                if isempty(find(lm_add_stock{add_i})) == 0 && isempty(find(lm_first_stock{add_i})) == 0
                    lm_add_init(add_i,:) = {first_gosa_obs{add_i,:} add_gosa_obs{add_i,:}};
                    lm_add_current(add_i,:) = {current_init_obs{add_i,:} current_add_obs{add_i,:}};
                    lm_add_size(add_i,:) = {lm_first_size{add_i,:} lm_current_size{add_i,:}};
                elseif isempty(find(lm_add_stock{add_i})) && isempty(find(lm_first_stock{add_i})) == 0
                    lm_add_init(add_i,:) = first_gosa_obs(add_i,:);
                    lm_add_current(add_i,:) = current_init_obs(add_i,:);
                    lm_add_size(add_i,:) = lm_first_size(add_i,:);
                elseif isempty(find(lm_add_stock{add_i}))  == 0&& isempty(find(lm_first_stock{add_i}))
                    lm_add_init(add_i,:) = add_gosa_obs(add_i,:);
                    lm_add_current(add_i,:) = current_add_obs(add_i,:);
                    lm_add_size(add_i,:) = lm_current_size(add_i,:);
                end
            end           
        end
                
        for obs_i = 1:length(glo_obs(1,:))
            r(obs_i)=en_plot_red(glo_obs(:,obs_i).',glo_rand_size(obs_i));
        end

        if isempty(find(sen_num == 0)) == 0
            [~,glo_range_obs]=sensor_judge(glo_obs,sen_num,glo_rand_size);
            [up_obs]=gosamodel(glo_range_obs,start,ang_wp);
            kill_sensor_range = plot_sensor_range(start,ang_wp);

            % 障害物の座標格納
            [~,cols] = size(up_obs);
            if cols < 2 && up_obs(2,1) == 0
                obs_list = [];
            else
                ob=ob_round(up_obs,rand_size);
                obs_list=ob.';
                
                %LMの変化量を図示
                L = lm_line(up_obs,gosa_obs);
                
                %ポテンシャル遷移で評価
                po_cdc_st=potential(up_obs,start,rand_size);
                po_cdc = [po_cdc po_cdc_st];
                sum_po_cdc=sum(po_cdc)/length(po_cdc);
                currentFile = sprintf('./potential/potential_cdc.mat');
                save(currentFile,'glo_obs','glo_gosa_obs','glo_rand_size','drive_cdc','po_cdc','sum_po_cdc');
                for j=1:length(up_obs(1,:))
                    b(j)=en_plot_orange(up_obs(:,j).',rand_size(j));
                end
            end
        end

        if strcmp(add_path.State,'finished') && flag_wp_continue == 0
            if isempty(add_path.Error)
                [add_path_stock,wp_add_stock] = fetchOutputs(add_path);
                wp_add(add_count,:) = {wp_add_stock(:,2:end)};
                disp('Path Add !!');
                plot(add_path_stock(1,:),add_path_stock(2,:),'Color',[1, 0, 0, 0.2],'LineWidth',3);
                hold on;
                cancel(add_path);
                flag_wp_continue = 1;
            end
        end
        
        if R == 0
            [wp,k,mat_er,plan_er,flag_stock]=compensation(up_obs,gosa_obs,lm_add_current,lm_add_init);
            
            if flag_stock == 1
                flag = 1;
            end
            
            if flag == 1
                ang = x(3);
                break;
            end
        end


        if isempty(find(select_add_num == 0)) == 0
            [add_range_size,add_range_obs]=sensor_judge(ground_add_obs,select_add_num,add_obs_rand_size);
            [up_add_obs]=gosamodel(add_range_obs,start,ang_wp);
            
            for obs_i = 1:length(up_add_obs(1,:))
                d_obs(obs_i) = circle_plot(up_add_obs(:,obs_i).',add_range_size(obs_i),mintblue);
            end
        end

        if add_flag == 1 && i > 1 
            delete(add_obs_plot);
            add_flag = 0;
        end
        
        
        if isempty(lm_add_stock) == 0
            for add_i = 1:add_count
                for obs_i = 1:length(lm_add_current{add_i}(1,:))
                    add_obs_plot_stock(obs_i) = en_plot_orange(lm_add_current{add_i}(:,obs_i).',lm_add_size{add_i}(obs_i));
                end
                add_obs_plot = [add_obs_plot add_obs_plot_stock];
            end
            add_flag = 1;
        end

        for j=wp_i:length(wp(1,:))
          wp_plt(j) = plot(wp(1,j),wp(2,j),'g:o','MarkerSize',10);
          hold on;
        end
        
        %{
        if rand_n == 9 && wp_i > 1
            disp('Add !!');
            flag = 1;
            flag_add = 1;
            ang = x(3);
            %プロットポイントコメントアウト部分
            delete(wp_plt);
            delete(L);
            delete(r);       
            break;
        end
        %}
        
        %ゴール判定
        if norm([x(1),x(2)]-goal) < Goal_tor
            disp('Arrive Goal!!');
            ang = x(3);
            %プロットポイントコメントアウト部分
            if length(wp(1,:)) == wp_i
                break;
            end
            delete(b);        
            break;
        end
        
        if i > 20 && norm([x(1) - result.x(end-20,1),x(2) - result.x(end-20,2)]) < 0.3
            flag = 0;
        end
        %{
        if  i > 150 && abs(result.x(length(result.x(:,1)),1) - result.x(length(result.x(:,1))-150,1)) < 1.0 && abs(result.x(length(result.x(:,1)),2) - result.x(length(result.x(:,1))-150,2)) < 1.0 
            disp('Skip Waypoint');
            %プロットポイントコメントアウト部分
               delete(b);
               delete(L);
               delete(r);
            break;
        end
        %}
        
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
         pause(0.01);
    end
    %プロットポイントコメントアウト部分
    delete(kill_sensor_range);
    delete(wp_plt);
    delete(L);
    delete(r);
    if i>1 & isempty(find(select_add_num == 0)) == 0
        delete(d_obs);
    end
    if add_flag == 1 && i > 1
        delete(add_obs_plot);
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
        [~,cols] = size(ob);
        if cols < 2 
            dist = 0;
        else
            dist=CalcDistEval(xt,ob,R);
        end
        
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
if x(3) > 2*pi
       ang = rem(x(3),2*pi);
else
    ang = x(3);
end
theta=toDegree(ang);%ロボットの方位
goalrad = atan2(goal(2)-x(2),goal(1)-x(1));
if goalrad < 0
    goalRad = 2*pi + goalrad;
else
    goalRad = goalrad;
end
goalTheta=toDegree(goalRad);%ゴールの方位

if goalTheta > 270 && theta + 360 - goalTheta < 180 ||  theta > 270 && goalTheta + 360 - theta < 90
    targetTheta = abs(theta + 360 - goalTheta);    
else
    targetTheta = abs(theta - goalTheta);    
end
heading=360-targetTheta;
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
 [X,Y] = pol2cart(x(3),dt); 
B = [X 0
    Y 0
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