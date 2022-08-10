currentFile = sprintf('path_interp_1_50.mat');
load(currentFile);
graph(glo_gosa_obs,glo_rand_size);

plot(drive_cdc(1,:),drive_cdc(2,:),'r:.');
hold on;

%初期宣言
Auto_path = drive_cdc;
obs = ob_round(glo_gosa_obs,glo_rand_size);
count = 1;
count_global =1;

WP_AUTO=DouglasPeucker(Auto_path,0.51);
while 1
Comand_path = [];
Comand_path_nosampling = [Auto_path(1,1) Auto_path(1,end)
                          Auto_path(2,1) Auto_path(2,end)];
L = linspace(Auto_path(2,1),Auto_path(2,end),length(Auto_path(1,:)));
Comand_path = interp1(Comand_path_nosampling(2,:),Comand_path_nosampling(1,:),L);
Comand_path = [Comand_path 
               L];
plot(Comand_path(1,:),Comand_path(2,:),'b:.');
hold on;
wp_auto = WP_AUTO;
wp_comand = DouglasPeucker(Comand_path,0.51);
plot(wp_auto(1,:),wp_auto(2,:),'g:o');
hold on;

if length(wp_auto(1,:)) == length(wp_comand(1,:))
    L = len(wp_auto(:,end).',drive_cdc(:,end).');
    if L < 0.5
        disp("start Navi !!");
        break;
    end
    WP_AUTO = [WP_AUTO(:,1:count_global) WP WP_AUTO(:,count_global+1:end)];
    count_global = count_global + count;
    Auto_path = drive_cdc;
    count = count_global + count;
    continue;
else
    if count == count_global
        WP = wp_auto(:,2:end);
    else
        WP =[WP wp_auto(:,2:end)];
    end
    p.start = [wp_auto(1,count),wp_auto(2,count)];
    p.goal = [wp_auto(1,count+1),wp_auto(2,count+1)];
    path = DynamicWindowApproach_global(p.start,p.goal,obs.');
    plot(path(1,:),path(2,:),'g:.');
    hold on;
    Auto_path = path;
    count = count + 1;
end
pause(0.0001);
end


function l=len(a,b)
l=norm(a-b);
end


function graph(glo_gosa_obs,size)
  for i=1:length(glo_gosa_obs(1,:))
    en_plot_blue(glo_gosa_obs(:,i).',size(i));
  end
    hold on;
    grid on;
    xlabel('x[m]')
    ylabel('y[m]')
    xlim([-20 20]);
    ylim([0 50]);
end

function en_plot_blue(glo_obs,size)
 [x,y]=circle(glo_obs(1,1),glo_obs(1,2),size);
 fill(x,y,'b');
 hold on;
end

function [r_x,r_y]=circle(x,y,r)
 t=linspace(0,2*pi,100);
 r_x=r*cos(t)+x;
 r_y=r*sin(t)+y;
end

function obstacle=ob_round(obs,r)
obstacle = [];
for i=1:length(r)
    [x,y]=circle(obs(1,i),obs(2,i),r(i));
    if i==1
        obstacle=[x;y];
    else
        obstacle=[obstacle(1,:) x;obstacle(2,:) y];
    end
end
end

