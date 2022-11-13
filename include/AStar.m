%% Astar path planning
function path = AStar(obstacle,start,goal)
% A*法によって最短経路を探索するプログラム
% 最短経路のパスの座標リストを返す

% start = [x
%          y];

% goal = [x
%         y];

% obstacle = [x ... x
%             y ... y];
p.start = start;
p.goal = goal;
oct_path=[];%パス
%計算中ノード情報格納用[x,y,cost,px(親ノード),py(親ノード)] startノードを格納する
open=[p.start(1,1) p.start(2,1) len(p.start.',p.goal.') p.start(1,1) p.start(2,1)];
close=[];%計算済みノード情報格納用

%隣接ノードへの移動モデル これを変えることでロボットの移動を指定できる
next=MotionModel();

findFlag=false;%ゴール発見フラグ

while ~findFlag
      %openにデータがない場合はパスが見つからなかった。
      
      if isempty(open)
          disp('No path to goal!!'); 
          path = [];
          return; 
      end

      %openなノードの中で最もコストが小さいものを選ぶ
      [Y,I] = sort(open(:,3));
      open=open(I,:);

      %ゴール判定
      if isSamePosi(open(1,1:2),p.goal.')
          disp('Find Goal!!');
          %ゴールのノードをCloseの先頭に移動
          close=[open(1,:);close];
          open(1,:)=[];
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
      %パス探索のステップ動画
      %animation(open,close,p,obstacle);
end

%最短パスの座標リストを取得
oct_path=GetPath(close,p.start.');
path=PathSmoothing(oct_path);
end

function result=isSamePosi(a,b)
%2x1のベクトルが同じかどうかを判断する関数
result=false;
com_x=abs(a(1)-b(1));
com_y=abs(a(2)-b(2));
l = sqrt(com_x^2+com_y^2);
if l < 0.5
    result=true;
end
end

function l=len(a,b)
l=norm(a-b);
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
next=[0.3 0.3 0.3
     0.3 0 0.3
      0 0.3 0.3
      -0.3 0 0.3
      0 -0.3 0.3
      -0.3 -0.3 0.3];
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
    flag=3;
    return;
end

function animation(open,close,p,obstacle)
% 探索の様子を逐次的に表示する関数

figure(1)
if length(obstacle)>=1
    plot(obstacle(:,1),obstacle(:,2),'om');hold on;
end
plot(p.start(1),p.start(2),'*r');hold on;
plot(p.goal(1),p.goal(2),'*b');hold on;
plot(open(:,1),open(:,2),'xr');hold on;
plot(close(:,1),close(:,2),'xk');hold on;
grid on;
pause(0.001);

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