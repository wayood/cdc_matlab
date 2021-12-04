load('potential.mat');%middleなし
load('potential_cdc.mat');%軌道補正込み
Potential_Field(glo_obs,glo_rand_size);
hold on;
for i=1:length(drive(1,:))
    z(i)=potential(glo_obs,drive(:,i),glo_rand_size);
    plot3(drive(1,i)+11,drive(2,i),z(i),'o','Color','b');
    hold on;
end
for j=1:length(drive_cdc(1,:))
    z_cdc(j)=potential(glo_obs,drive_cdc(:,j),glo_rand_size);
    plot3(drive_cdc(1,j)+11,drive_cdc(2,j),z(j),'o','Color','r');
    hold on;
end
pause();
hold off;
i=1:length(drive(1,:));
plot(i,z(i));
hold on;
i=1:length(drive_cdc(1,:));
plot(i,z_cdc(i));
grid on;
xlabel('x[m]');
ylabel('potential');

function po=potential(obs,move,size)
    po=0;
    for i=1:obs(1,:)+1
     l=norm(obs(:,i).'-move.');
     if l < size(i)
       p=(3-l.^2/size(i).^2)/2;
     else
       p=size(i)/l;
     end
       po=po+p;
    end
end