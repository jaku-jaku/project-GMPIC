%% Utilities %%
classdef utils
    methods(Static)
        %% Plotting %%
        function plot_link( ...
            curret_config_sframe_SE3_base, ...
            link_relative_displacement_R3, ...
            link_base_frame_name, ...
            link_color, style)
            %{
                Plot the link per given the following:
                base |/_ x---->O tip |/_
            %}
            origin = [0;0;0;1]; % w.r.t the current link
            end_pos = [link_relative_displacement_R3;1];
            % Homogeneous form
            points_on_link = [origin end_pos];
            
            % transform link wrt spatial base frame
            g_tip = curret_config_sframe_SE3_base * points_on_link;
            % Finally, we plot the points
            plot3(g_tip(1,:),g_tip(2,:),g_tip(3,:),'-','Color',link_color,'LineWidth',4)
            hold on
            % end
            % Plot the co-ordinate frame of the link at the base:
            utils.plot_axis(curret_config_sframe_SE3_base,0.1,style,link_base_frame_name)
        end
        function plot_axis_screw_motion( ...
            screw_motion_SE3, ...
            spatial_frame_SE3, ...
            scale, style, ...
            axis_name )
            % Plot the co-ordinate frame of the link
            origin = [0;0;0;1]; % w.r.t the current link
            axis_coordinate_frame = [scale 0 0;...
                0 scale 0;...
                0 0 scale;...
                1 1 1];
            world_axis_coordinate_frame = screw_motion_SE3*spatial_frame_SE3*axis_coordinate_frame;
            transformed_origin = screw_motion_SE3*spatial_frame_SE3*origin;
            text(transformed_origin(1,1),transformed_origin(2,1),transformed_origin(3,1),axis_name)
            x_axis = [transformed_origin world_axis_coordinate_frame(:,1)];
            y_axis = [transformed_origin world_axis_coordinate_frame(:,2)];
            z_axis = [transformed_origin world_axis_coordinate_frame(:,3)];
            % Plot x-axis
            plot3(x_axis(1,:),x_axis(2,:),x_axis(3,:),sprintf('%sr', style),'LineWidth',2)
            hold on
            % Plot y-axis
            plot3(y_axis(1,:),y_axis(2,:),y_axis(3,:),sprintf('%sg', style),'LineWidth',2)
            hold on
            % Plot z-axis
            plot3(z_axis(1,:),z_axis(2,:),z_axis(3,:),sprintf('%sb', style),'LineWidth',2)
            hold on
            % Plot the coordinate frame of the tail:
        end
        function plot_axis( ...
            spatial_frame_SE3, ...
            scale, style, ...
            axis_name )
            % Plot the co-ordinate frame of the link
            origin = [0;0;0;1]; % w.r.t the current link
            axis_coordinate_frame = [scale 0 0;...
                0 scale 0;...
                0 0 scale;...
                1 1 1];
            world_axis_coordinate_frame = spatial_frame_SE3*axis_coordinate_frame;
            transformed_origin = spatial_frame_SE3*origin;
            text(transformed_origin(1,1),transformed_origin(2,1),transformed_origin(3,1),axis_name)
            x_axis = [transformed_origin world_axis_coordinate_frame(:,1)];
            y_axis = [transformed_origin world_axis_coordinate_frame(:,2)];
            z_axis = [transformed_origin world_axis_coordinate_frame(:,3)];
            % Plot x-axis
            plot3(x_axis(1,:),x_axis(2,:),x_axis(3,:),sprintf('%sr', style),'LineWidth',2)
            hold on
            % Plot y-axis
            plot3(y_axis(1,:),y_axis(2,:),y_axis(3,:),sprintf('%sg', style),'LineWidth',2)
            hold on
            % Plot z-axis
            plot3(z_axis(1,:),z_axis(2,:),z_axis(3,:),sprintf('%sb', style),'LineWidth',2)
            hold on
            % Plot the coordinate frame of the tail:
        end
        %% Plotting Summit %%
        function transformed_origin = plot_Summit(screw_rotation,loc_of_link_frame,link_name,axis_ref)
            %UNTITLED4 PLot the summit
            % screw_rotation - This gives us the total solid-motion undergone by the
            % link
            % loc_of_link_frame - This is the initial location of the origin and
            % orientation of the link. The above screw motion acts on this link
            % loc_of_next_link_frame - This gives us just the location of where the
            % next successive link starts from the origin of the current link, w.r.t
            % the current link
            origin = [0;0;0;1]; % w.r.t the current link
            if(link_name=="S")
                length = 0.62;
                width = 0.43;
                height = 0.56;
                points = [length/2 width/2 height 1;...
                          length/2 -width/2 height 1;...
                          -length/2 -width/2 height 1;...
                          -length/2 width/2 height 1];
                color = [158/255 158/255 158/255];
            else
                length = 0.84;
                width = 0.48;
                height = 0.54;
                points = [ width/2 -height length/2 1;...
                           -width/2 -height length/2 1;...
                           -width/2 -height -length/2 1;...
                           width/2 -height -length/2 1];
                color = [0 0.4470 0.7410];
            end
            % Homogeneous form
            points_on_link = points';
            % Now, we transform the points from the link-frame to the world frame
            % before rigid body motion
            before_rigid_body_motion = loc_of_link_frame*points_on_link;
            % Now we transform the body from its initial position to its final position
            after_rigid_body_motion = screw_rotation*before_rigid_body_motion;
            % Finally, we plot the points
            link = after_rigid_body_motion(:,1:2);
            plot3(link(1,:),link(2,:),link(3,:),'-','Color',color,'LineWidth',6)
            % if(~isempty(axis_ref))
            %     hold(axis_ref,'on')
            % else
            %     hold on
            % end
            hold on
            link = after_rigid_body_motion(:,2:3);
            plot3(link(1,:),link(2,:),link(3,:),'-','Color',color,'LineWidth',6)
            % if(~isempty(axis_ref))
            %     hold(axis_ref,'on')
            % else
            %     hold on
            % end
            hold on
            link = after_rigid_body_motion(:,3:4);
            plot3(link(1,:),link(2,:),link(3,:),'-','Color',color,'LineWidth',6)
            % if(~isempty(axis_ref))
            %     hold(axis_ref,'on')
            % else
            %     hold on
            % end
            hold on
            link = [after_rigid_body_motion(:,4) after_rigid_body_motion(:,1)];
            plot3(link(1,:),link(2,:),link(3,:),'-','Color',color,'LineWidth',6)
            % if(~isempty(axis_ref))
            %     hold(axis_ref,'on')
            % else
            %     hold on
            % end
            hold on
            % if(link_name=="C")
            %     3
            % end
            % Plot the co-ordinate frame of the link at the base:
            utils.plot_axis_screw_motion(screw_rotation,loc_of_link_frame,0.1,'-',link_name)
        end
    end
end