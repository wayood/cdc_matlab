function [wp_add] = voronoi_waypoint_generation(start,goal,obstacle)

    % start and goal is [x
    %                    y];
    
    % obstacle [x … x
    %           y … y];
    
    obstacle = [start obstacle goal];
    P = [obstacle(1,:).' obstacle(2,:).'];
    DT = delaunayTriangulation(P);
    [wp_canditate,edge] = voronoiDiagram(DT);
    [vx,vy] = voronoi(obstacle(1,:),obstacle(2,:));
    voronoi(obstacle(1,:),obstacle(2,:));
    hold on;

    % スタートとゴール地点での周辺のボロノイ点を代入
    
    start_find_vx = wp_canditate(edge{1},1);
    end_find_vx = wp_canditate(edge{length(obstacle(1,:))},1);
    
    vx_num = ones(4,length(vx(1,:)));
    vy_num = ones(4,length(vx(1,:)));
    find_cross_x = zeros(4,1);
    find_cross_y = zeros(4,1);
    weights = [];
    num_count = 1;
    
    for i = 1:length(vx(1,:))
        if find(vx(1,i) == vx_num & vy(1,i) == vy_num)
            [cols,rows] = find(vx(1,i) == vx_num & vy(1,i) == vy_num);
            if cols(1,1) > 2
                cols(1,1) = cols(1,1) - 2;
            end
            vx_num(1,i) = vx_num(cols(1,1),rows(1,1));
            vx_num(3,i) = vx(1,i);
            vy_num(1,i) = vy_num(cols(1,1),rows(1,1));
            vy_num(3,i) = vy(1,i);
        else
            vx_num(1,i) = num_count;
            vx_num(3,i) = vx(1,i);
            vy_num(1,i) = num_count;
            vy_num(3,i) = vy(1,i);
            num_count = num_count + 1;
        end

        if find(vx(2,i) == vx_num & vy(2,i) == vy_num)
            [cols,rows] = find(vx(2,i) == vx_num & vy(2,i) == vy_num);
            if cols(1,1) > 2
                cols(1,1) = cols(1,1) -2;
            end
            vx_num(2,i) = vx_num(cols(1,1),rows(1,1));
            vx_num(4,i) = vx(2,i);
            vy_num(2,i) = vy_num(cols(1,1),rows(1,1));
            vy_num(4,i) = vy(2,i);
        else
            vx_num(2,i) = num_count;
            vx_num(4,i) = vx(2,i);
            vy_num(2,i) = num_count;
            vy_num(4,i) = vy(2,i);
            num_count = num_count + 1;
        end
        
        weights = [weights norm([vx_num(3,i) vy_num(3,i)]-[vx_num(4,i) vy_num(4,i)])];
        
        
        [pol_x,pol_y] = polyxpoly([start(1,1) goal(1,1)],[start(2,1) goal(2,1)],[vx_num(3,i) vx_num(4,i)],[vy_num(3,i) vy_num(4,i)]);
        
        if isempty(pol_x) == 0 && isempty(pol_y) == 0
            stock = [num_count num_count
                     vx_num(1,i) vx_num(2,i)
                     pol_x pol_x
                     vx_num(3,i) vx_num(4,i)];
                 
            find_cross_x = [find_cross_x stock];
            
            stock = [num_count num_count
                     vy_num(1,i) vy_num(2,i)
                     pol_y pol_y
                     vy_num(3,i) vy_num(4,i)];
                 
            find_cross_y = [find_cross_y stock];
            
            num_count = num_count + 1;
        end
    end
           
    
    for i = 1:length(find_cross_x(1,:))-1
        if size(start_find_vx) < 3
            if find(start_find_vx == find_cross_x(4,i) | start_find_vx == find_cross_x(4,i+1))
                vx_num = [vx_num find_cross_x(:,i:i+1)];
                vy_num = [vy_num find_cross_y(:,i:i+1)];
                vx = [vx find_cross_x(3:4,i:i+1)];
                vy = [vy find_cross_y(3:4,i:i+1)];
                weights = [weights norm([find_cross_x(3,i) find_cross_y(3,i)]-[find_cross_x(4,i) find_cross_y(4,i)]) norm([find_cross_x(3,i+1) find_cross_y(3,i+1)]-[find_cross_x(4,i+1) find_cross_y(4,i+1)])];
                start_node_num = find_cross_x(1,i);
                plot(find_cross_x(3,i),find_cross_y(3,i),"g:*");
                hold on;
            end
        end
        
        if find(start_find_vx == find_cross_x(4,i))
            if find(start_find_vx == find_cross_x(4,i+1))
                vx_num = [vx_num find_cross_x(:,i:i+1)];
                vy_num = [vy_num find_cross_y(:,i:i+1)];
                vx = [vx find_cross_x(3:4,i:i+1)];
                vy = [vy find_cross_y(3:4,i:i+1)];
                weights = [weights norm([find_cross_x(3,i) find_cross_y(3,i)]-[find_cross_x(4,i) find_cross_y(4,i)]) norm([find_cross_x(3,i+1) find_cross_y(3,i+1)]-[find_cross_x(4,i+1) find_cross_y(4,i+1)])];
                start_node_num = find_cross_x(1,i);
                plot(find_cross_x(3,i),find_cross_y(3,i),"g:*");
                hold on;
            end
        end
        
        if find(end_find_vx == find_cross_x(4,i))
            if find(end_find_vx == find_cross_x(4,i+1))
                vx_num = [vx_num find_cross_x(:,i:i+1)];
                vy_num = [vy_num find_cross_y(:,i:i+1)];
                vx = [vx find_cross_x(3:4,i:i+1)];
                vy = [vy find_cross_y(3:4,i:i+1)];
                weights = [weights norm([find_cross_x(3,i) find_cross_y(3,i)]-[find_cross_x(4,i) find_cross_y(4,i)]) norm([find_cross_x(3,i+1) find_cross_y(3,i+1)]-[find_cross_x(4,i+1) find_cross_y(4,i+1)])];
                goal_node_num = find_cross_x(1,i);
                plot(find_cross_x(3,i),find_cross_y(3,i),"k:*");
                hold on;
            end
        end
        
    end
        
        
    s = vx_num(1,:);
    t = vx_num(2,:);
    G = graph(s,t,weights);

    [path,~] = shortestpath(G,start_node_num,goal_node_num);
    wp_add = [];
    
    for i = 1:length(path)
        [cols,rows] = find(path(1,i) == vx_num);
        wp = [vx(cols(1,1),rows(1,1))
              vy(cols(1,1),rows(1,1))];
        wp_add = [wp_add wp];
    end
    wp_add = [wp_add goal];
    
end


