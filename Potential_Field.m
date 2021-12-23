function []=Potential_Field(obs,R,x,y)

x_p=-x/2:x/2;
y_p=0:y;
for i=1:x
    for j=1:y
        s=[x_p(i);y_p(j)];
        z(j,i)=potential(obs,s,R);
    end
end
surf(z);
xlabel('x[m]');
ylabel('y[m]');
zlabel('potential');

end

function po=potential(obs,move,size)
  po=0;
    for i=1:obs(1,:)+1
     l=norm(obs(:,i).'- move.');
     if l < size(i)
       p=(3-l.^2/size(i).^2)/2;
     else
       p=size(i)/l;
     end
       po=po+p;
    end
end

%一雇用のポテンシャル場
function po=potential_1(obs,move,size)
    po=0;
    l=norm(obs.'-move.');
    if l < size
       p=(3-l.^2/size.^2)/2;
    else
       p=size/l;
    end
    po=po+p;    
end