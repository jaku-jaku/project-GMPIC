%% Utilities %%
classdef utils
    methods(Static)
        %% Plotting %%
        function plot_link(screw_rotation,loc_of_link_frame,loc_of_next_link_frame,link_name,link_color,axis_ref)
            %UNTITLED Given the location and orientation of the link and the link
            %dimension, plot it
            %   Detailed explanation goes here
            % screw_rotation - This gives us the total solid-motion undergone by the
            % link
            % loc_of_link_frame - This is the initial location of the origin and
            % orientation of the link. The above screw motion acts on this link
            % loc_of_next_link_frame - This gives us just the location of where the
            % next successive link starts from the origin of the current link, w.r.t
            % the current link
            % TODO: revise this function to have link tip frame plotted instead
            origin = [0;0;0;1]; % w.r.t the current link
            end_pos = [loc_of_next_link_frame;1];
            if(link_name=="W3")
                % Homogeneous form
                points_on_link = [origin [0;0;loc_of_next_link_frame(3,1);1] end_pos];
            elseif(link_name=="W4")
                % Homogeneous form
                points_on_link = [origin [loc_of_next_link_frame(1,1);0;0;1] end_pos];
            else
                % Homogeneous form
                points_on_link = [origin end_pos];
            end
            % Now, we transform the points from the link-frame to the world frame
            % before rigid body motion
            before_rigid_body_motion = loc_of_link_frame*points_on_link;
            % Now we transform the body from its initial position to its final position
            after_rigid_body_motion = screw_rotation*before_rigid_body_motion;
            % Finally, we plot the points
            plot3(after_rigid_body_motion(1,:),after_rigid_body_motion(2,:),after_rigid_body_motion(3,:),'-','Color',link_color,'LineWidth',4)
            % if(link_name=='S2')
            %     after_rigid_body_motion
            % end
            % if(~isempty(axis_ref))
            %     hold(axis_ref,'on')
            % else
            %     hold on
            % end
            hold on
            % end
            % Plot the co-ordinate frame of the link
            body_link_coordinate_frame = [0.1 0 0;...
                                          0 0.1 0;...
                                          0 0 0.1;...
                                          1 1 1];
            world_link_coordinate_frame = screw_rotation*loc_of_link_frame*body_link_coordinate_frame;
            transformed_origin = screw_rotation*loc_of_link_frame*origin;
            text(transformed_origin(1,1),transformed_origin(2,1),transformed_origin(3,1),link_name)
            x_axis = [transformed_origin world_link_coordinate_frame(:,1)];
            y_axis = [transformed_origin world_link_coordinate_frame(:,2)];
            z_axis = [transformed_origin world_link_coordinate_frame(:,3)];
            % if(link_name=='Summit')
            %     world_link_coordinate_frame
            % end
            % Plot x-axis
            plot3(x_axis(1,:),x_axis(2,:),x_axis(3,:),'-r','LineWidth',2)
            % if(~isempty(axis_ref))
            %     hold(axis_ref,'on')
            % else
            %     hold on
            % end
            hold on
            % Plot y-axis
            plot3(y_axis(1,:),y_axis(2,:),y_axis(3,:),'-g','LineWidth',2)
            % if(~isempty(axis_ref))
            %     hold(axis_ref,'on')
            % else
            %     hold on
            % end
            hold on
            % Plot z-axis
            plot3(z_axis(1,:),z_axis(2,:),z_axis(3,:),'-b','LineWidth',2)
            % if(~isempty(axis_ref))
            %     hold(axis_ref,'on')
            % else
            %     hold on
            % end
            hold on
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
            % Plot the co-ordinate frame of the link
            body_link_coordinate_frame = [0.1 0 0;...
                                          0 0.1 0;...
                                          0 0 0.1;...
                                          1 1 1];
            world_link_coordinate_frame = screw_rotation*loc_of_link_frame*body_link_coordinate_frame;
            transformed_origin = screw_rotation*loc_of_link_frame*origin;
            text(transformed_origin(1,1),transformed_origin(2,1),transformed_origin(3,1),link_name)
            x_axis = [transformed_origin world_link_coordinate_frame(:,1)];
            y_axis = [transformed_origin world_link_coordinate_frame(:,2)];
            z_axis = [transformed_origin world_link_coordinate_frame(:,3)];
            
            % if(link_name=='Summit')
            %     world_link_coordinate_frame
            % end
            % Plot x-axis
            plot3(x_axis(1,:),x_axis(2,:),x_axis(3,:),'-r','LineWidth',2)
            % if(~isempty(axis_ref))
            %     hold(axis_ref,'on')
            % else
            %     hold on
            % end
            hold on
            % Plot y-axis
            plot3(y_axis(1,:),y_axis(2,:),y_axis(3,:),'-g','LineWidth',2)
            % if(~isempty(axis_ref))
            %     hold(axis_ref,'on')
            % else
            %     hold on
            % end
            hold on
            % Plot z-axis
            plot3(z_axis(1,:),z_axis(2,:),z_axis(3,:),'-b','LineWidth',2)
            % if(~isempty(axis_ref))
            %     hold(axis_ref,'on')
            % else
            %     hold on
            % end
            hold on
        end
    end
end