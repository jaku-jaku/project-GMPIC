classdef common
    properties(Constant)
        %% [ WAM constant ]:
        % - numerical:
        L_SHOULDER    =0.346;
        L_ARM         =0.55;
        L_ELBOW       =0.3;
        L_WRIST       =0.06;
        W_ELBOW_JNT   =0.045;
        % - symbolic:
        % L_SHOULDER    = sym("l0");
        % L_ARM         = sym("l1");
        % L_ELBOW       = sym("l2");
        % W_ELBOW_JNT   = sym("w1");
        %% [ SUMMIT Model ]:
        SUMMIT_INIT_POSE = eye(4);
        G_SUMMIT_WAM = [ 0.14 0 0.405 ]; 
        %% [ WAM Model ]:
        RELATIVE_Links = [
            [ 0 0 common.L_SHOULDER ];                % 0-1
            [ 0 0 0 ];                                % 1-2
            [ 0 0 0 ];                                % 2-3
            [ common.W_ELBOW_JNT 0 common.L_ARM ];    % 3-4
            [ -common.W_ELBOW_JNT -common.L_ELBOW 0 ];% 4-5
            [ 0 0 0 ];                                % 5-6
            [ 0 0 0 ];                                % 6-7
            [ 0 0 0.2]; % 7-8
        ]; % N+1 links
        RELATIVE_Axes = [
            [ 0 0 1 ]; % 0-1
            [ 1 0 0 ]; % 1-2
            [ 1 0 0 ]; % 2-3
            [ 1 0 0 ]; % 3-4
            [ 1 0 0 ]; % 4-5
            [ 1 0 0 ]; % 5-6
            [ 1 0 0 ]; % 6-7
        ]; % N jnts
        RELATIVE_Axes_Rotations = [
            0;           % 0-1
            -pi/2;       % 1-2
            pi/2;        % 2-3
            -pi/2;       % 3-4
            pi/2;        % 4-5
            -pi/2;       % 5-6
            pi/2;        % 6-7
        ]; % N jnts
        % - configs:
        % INIT_JOINT_HOME = [pi;pi/2;0;0;0;0;0];
    end
    methods(Static)
        % TBD
    end
end