function [wp_add] = voronoi_waypoint_generation(start,goal,obstacle)

    % start and goal is 2cols,1rows;
    % obstacle [x … x
    %           y … y];
    
    obstacle = [start obstacle goal];
    P = [obstacle(1,:).' obstacle(2,:).'];
    DT = delaunayTriangulation(P);
    [wp_canditate,edge] = voronoiDiagram(DT);
    [vx,vy] = voronoi(obstacle(1,:),obstacle(2,:));

    x_start = wp_canditate(edge{1},1);
    y_start = wp_canditate(edge{1},2);
    x_end = wp_canditate(edge{length(obstacle(1,:))},1);
    y_end = wp_canditate(edge{length(obstacle(1,:))},2);

    start = length(vx) + 1;
    goal = start + length(x_start);

    for i = 1:length(x_start)
        x = [obstacle(1,1) 
             x_start(i,1)];
        vx = [vx x];
        y = [obstacle(2,1) 
             y_start(i,1)];
        vy = [vy y];
    end

    for i = 1:length(x_end)
        x = [obstacle(1,end) 
             x_end(i,1)];
        vx = [vx x];
        y = [obstacle(2,end) 
             y_end(i,1)];
        vy = [vy y];
    end

    vx_num = ones(4,length(vx(1,:)));
    vy_num = ones(4,length(vx(1,:)));
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

    end

    start_num = vx_num(1,start);
    goal_num = vx_num(1,goal);
    s = vx_num(1,:);
    t = vx_num(2,:);
    G = graph(s,t,weights);

    [path,~] = shortestpath(G,start_num,goal_num);
    wp_add = [];

    for i = 1:length(path)
        [cols,rows] = find(path(1,i) == vx_num);
        wp = [vx(cols(1,1),rows(1,1))
              vy(cols(1,1),rows(1,1))];
        wp_add = [wp_add wp];
    end
end


