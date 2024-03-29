%% 軌道補正アルゴリズム
global cdc_length;
cdc_length=0.5;%m単位で補正ポイントを設定

%% パラメータ設定
p.start=[0;0];
p.goal=[0;20];
for i=1:50
    ran_x=-10+20*rand;
    ran_y=20*rand;
    glo_rand_size(i)=0.3+0.5*rand;
    glo_obs(1,i)=ran_x;
    glo_obs(2,i)=ran_y;
end

[glo_gosa_obs]=gosamodel(glo_obs,p);%経路生成時のLM座標
obs=ob_round(glo_gosa_obs,glo_rand_size);
path=AStar(obs,p);
graph(glo_gosa_obs,p,glo_rand_size);
path=PathSmoothing(path);
if length(path)>=1
    plot(path(:,1),path(:,2),'-r','MarkerSize',3);
    hold on;
end
save('path_interp.mat','path','glo_obs','glo_gosa_obs','glo_rand_size','obs');
hold on;
for io=1:length(glo_obs(1,:))
  en_plot_red(glo_obs(:,io).',glo_rand_size(io));
  hold on;
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
    xlim([-10 10]);
    ylim([0 20]);
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
%% Astar path planning
function path=AStar(obstacle,p)
% A*法によって最短経路を探索するプログラム
% 最短経路のパスの座標リストを返す

path=[];%パス
%計算中ノード情報格納用[x,y,cost,px(親ノード),py(親ノード)] startノードを格納する
open=[p.start(1,1) p.start(2,1) len(p.start.',p.goal.') p.start(1,1) p.start(2,1)];
close=[];%計算済みノード情報格納用

%隣接ノードへの移動モデル これを変えることでロボットの移動を指定できる
next=MotionModel();

findFlag=false;%ゴール発見フラグ

while ~findFlag
      %openにデータがない場合はパスが見つからなかった。
      if isempty(open(:,1))
          disp('No path to goal!!'); 
          return; 
      end
      %openなノードの中で最もコストが小さいものを選ぶ
      [Y,I] = sort(open(:,3));
      open=open(I,:);
      
      %ゴール判定
      if isSamePosi(open(1,1:2),p.goal.')
          disp('Find Goal!!');
          %ゴールのノードをCloseの先頭に移動
          close=[open(1,:);close];open(1,:)=[];
          findFlag=true;
          break;
      end
      
      for in=1:length(next(:,1))
          %隣接ノードの位置とコストの計算
          m=[open(1,1)+next(in,1) open(1,2)+next(in,2) open(1,3)];
          m(3)=m(3)+1+len(m(1:2),p.goal.')-len(open(1,1:2),p.goal.');%コストの計算
          
          %隣接ノードが障害物だったら次のノードを探す
          if isObstacle(m,obstacle) 
              continue; 
          end
          
          %openとcloseのリストの中にmが含まれるかを探索
          [flag, targetInd]=FindList(m,open,close);

          if flag==1 %openリストにある場合
              if m(3)<open(targetInd,3)
                  open(targetInd,3)=m(3);
                  open(targetInd,4)=open(1,1);
                  open(targetInd,5)=open(1,2);
              end
          elseif flag==2 %closeリストにある場合
              if m(3)<close(targetInd,3)
                  %親ノードの更新
                  close(targetInd,4)=open(1,1);
                  close(targetInd,5)=open(1,2);
                  open=[open; close(targetInd,:)];
                  close(targetInd,:)=[];%Openリストに移動
              end
          else %どちらにも無い場合
              %openリストに親ノードのインデックスと共に追加
              open=[open;[m open(1,1) open(1,2)]];
          end
      end

      %隣接ノード計算したopenノードはcloseノードへ移動
      if findFlag==false
          close=[close; open(1,:)];
          open(1,:)=[];
      end
      
end

%最短パスの座標リストを取得
path=GetPath(close,p.start.');

end

function result=isSamePosi(a,b)
%2x1のベクトルが同じかどうかを判断する関数
result=false;
com_x=abs(a(1)-b(1));
com_y=abs(a(2)-b(2));
if com_x<=0.2 && com_y<=0.2
    result=true;
end
end


function flag=isObstacle(m,obstacle)

for io=1:length(obstacle(1,:))
    if isSamePosi(obstacle(:,io).',m(1:2))
        flag=true;
        return;
    end
end
flag=false;%障害物ではない
end

function next=MotionModel()
%隣接ノードへの移動モデル これを変えることでロボットの移動を指定できる
next=[0.3 0.3
     0.3 0
      0 0.3
      -0.3 0
      0 -0.3
      -0.3 -0.3];
end

function path=GetPath(close,start)
%スタートからゴールまでの座標リストを取得する関数
ind=1;%goalはcloseリストの先頭に入っている
path=[];
while 1
    %座標をリストに登録
    path=[path; close(ind,1:2)];
    
    %スタート地点まで到達したか判断
    if isSamePosi(close(ind,1:2),start)   
        break;
    end
    
    %closeリストの中で親ノードを探す
    for io=1:length(close(:,1))
        if isSamePosi(close(io,1:2),close(ind,4:5))
            ind=io;
            break;
        end
    end
end

end

function [flag, targetInd]=FindList(m,open,close)
    targetInd=0;
    %openリストにあるか?
    if ~isempty(open)
        for io=1:length(open(:,1))
            if isSamePosi(open(io,:),m(1:2))
                flag=1;
                targetInd=io;
                return;
            end
        end
    end
    %closeリストにあるか?
    if ~isempty(close)
        for ic=1:length(close(:,1))
            if isSamePosi(close(ic,:),m(1:2))
                flag=2;
                targetInd=ic;
                return;
            end
        end
    end
    %どちらにも無かった
    flag=3;return;
end

%% path 平滑化

function optPath=PathSmoothing(path)
optPath=path;%元のパスをコピー

%平準化パラメータ
alpha=0.5;
beta=0.2;

torelance=0.00001;%パスの変化量の閾値(変化量がこの値以下の時平滑化を終了)
change=torelance;%パスの位置の変化量
while change>=torelance 
    change=0;%初期化
    for ip=2:(length(path(:,1))-1) %始点と終点は固定
        prePath=optPath(ip,:);%変化量計測用
        optPath(ip,:)=optPath(ip,:)-alpha*(optPath(ip,:)-path(ip,:));
        optPath(ip,:)=optPath(ip,:)-beta*(2*optPath(ip,:)-optPath(ip-1,:)-optPath(ip+1,:));
        change=change+norm(optPath(ip,:)-prePath);
    end
end

end
