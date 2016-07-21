function [axis_limits] = DisplayPoints3D(Model, Scene, sampling, axis_limits, transformed_ctrl_pts, match_result)
%%=====================================================================
%% $RCSfile: DisplayPoints3D.m,v $
%% $Author: bing.jian $
%% $Date: 2008-11-13 16:34:29 -0500 (Thu, 13 Nov 2008) $
%% $Revision: 109 $
%%=====================================================================

set(gca,'FontSize',16,'FontName','Times','FontWeight','bold');

markerArray = '+o*.xsd^v><ph';

for i = 1:size(Model,1)
    %plot3(Model(i,1),Model(i,2),Model(i,3),'r+', 'MarkerSize', 5, 'LineWidth',1.5);
    plot3(Model(i,1),Model(i,2),Model(i,3),strcat('r', markerArray(mod(i-1,13)+1)), 'MarkerSize', 15, 'LineWidth',1.5);
    hold on;
end
for i = 1:size(Scene,1)
    %plot3(Scene(i,1),Scene(i,2),Scene(i,3),'bo', 'MarkerSize', 5, 'LineWidth',1.5);
    plot3(Scene(i,1),Scene(i,2),Scene(i,3),strcat('b', markerArray(mod(i-1,13)+1)), 'MarkerSize', 15, 'LineWidth',1.5);
    hold on;
end
if nargin > 4
    %control points or deformed control points are entered, display each
    %transformation step
%     for level = 2:size(transformed_ctrl_pts,3)
%         for i = 1:size(transformed_ctrl_pts,1)
%             %plot3(ctrl_vector(i,1),ctrl_vector(i,2),ctrl_vector(i,3),strcat('y.'), 'MarkerSize', 10, 'LineWidth',1);
%             quiver3(transformed_ctrl_pts(i,1,level-1),transformed_ctrl_pts(i,2,level-1),transformed_ctrl_pts(i,3,level-1), transformed_ctrl_pts(i,1,level) - transformed_ctrl_pts(i,1,level-1), transformed_ctrl_pts(i,2,level) - transformed_ctrl_pts(i,2,level-1), transformed_ctrl_pts(i,3,level) - transformed_ctrl_pts(i,3,level-1),'r','AutoScale','off');
%             hold on;
%         end
%     end

    %display the single transformation step
    level = size(transformed_ctrl_pts,3);
    for i = 1:size(transformed_ctrl_pts,1)
        %plot3(ctrl_vector(i,1),ctrl_vector(i,2),ctrl_vector(i,3),strcat('y.'), 'MarkerSize', 10, 'LineWidth',1);
        quiver3(transformed_ctrl_pts(i,1,1),transformed_ctrl_pts(i,2,1),transformed_ctrl_pts(i,3,1), transformed_ctrl_pts(i,1,level) - transformed_ctrl_pts(i,1,1), transformed_ctrl_pts(i,2,level) - transformed_ctrl_pts(i,2,1), transformed_ctrl_pts(i,3,level) - transformed_ctrl_pts(i,3,1),'r','AutoScale','off');
        hold on;
    end
    
    for i = 1:size(match_result,1)
        if match_result(i,1) ~= 0
            %there is a match, draw a yellow line to connect the match
            X = [Scene(i,1), Model(match_result(i,1),1)];
            Y = [Scene(i,2), Model(match_result(i,1),2)];
            Z = [Scene(i,3), Model(match_result(i,1),3)];
            plot3(X, Y, Z , 'y', 'LineWidth', 1.5);
            hold on;
        end
    end
end
hold off;

axis equal;

if (nargin<3)
    sampling = 0;
end

m = size(Model,1);
if (sampling>0)
    for i=1:sampling:m
        text(Model(i,1), Model(i,2), Model(i,3), [' \leftarrow',sprintf('%d',i)]);
    end
end

m = size(Scene,1);
if (sampling>0)
    for i=1:sampling:m
        text(Scene(i,1), Scene(i,2), Scene(i,3), [' \leftarrow',sprintf('%d',i)]);
    end
end

if (nargin<4)
    axis_limits = determine_border(Model, Scene);
end

xlim(axis_limits(1,:));
ylim(axis_limits(2,:));
zlim(axis_limits(3,:));

%pbaspect([1,1,1]);
