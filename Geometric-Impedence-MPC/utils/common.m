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
        dP_SUMMIT_WAM = [ 0.14 0 0.405 ]; 
        % - configs:
        INIT_JOINT_HOME = [0 0 0 0 0 0 0];
        %% [ WAM Model ]:
        % - [ Dynamics from DataSheet ]:
        WAM = [...
            struct( ... % Link 1: base: 0--->1 :J1
                'name',         'base_B3350', ...
                'base_frame',   'F0', ...
                'mass',         9.97059584, ...
                'mass_center',  [-0.02017671, -0.26604706, -0.14071720], ...
                'MoI_at_mass',  [   0.10916849, 0.00640270, 0.02557874;
                                    0.00640270, 0.18294303, 0.00161433;
                                    0.02557874, 0.00161433, 0.11760385], ...
                'tip_frame',    'F1', ...
                'tip_frame_transformation', struct( ...
                    'q', [ 0 0 common.L_SHOULDER ], ...
                    'w', [ 0 0 1 ], ...
                    't', 0 ...
                ), ...
                'tip_frame_revolute_joint', struct( ...
                    'name', 'J1', ...
                    'axis', [ 0 0 1 ], ... % axis of roation at z-axis
                    'limits', [ -2.6, 2.6 ] ...
                ) ...
            ), ...
            struct( ... % 1-->2
                'name',                 'shoulder_B3351', ...
                'base_frame',           'F1', ...
                'mass',                 10.76768767, ...
                'mass_center',          [-0.00443422, 0.12189039, -0.00066489], ...
                'MoI_at_mass',          [   0.13488033, -0.00213041, -0.00012485; 
                                            -0.00213041,  0.11328369,  0.00068555;
                                            -0.00012485,  0.00068555,  0.09046330 ], ...
                'tip_frame',            'F2', ...
                'tip_frame_transformation', struct( ...
                    'q', [ 0 0 0 ], ...
                    'w', [ 1 0 0 ], ...
                    't', -pi/2 ...
                ), ...
                "tip_frame_revolute_joint", struct( ...
                    'name', 'J2', ...
                    'axis', [ 0 0 1 ], ... % axis of roation at z-axis
                    'limits', [ -2.0, 2.0 ] ...
                ) ...
            ), ...
            struct( ... 
                'name',                 'arm_base_B3352', ...
                'base_frame',           'F2', ...
                'mass',                 3.87493756, ...
                'mass_center',          [-0.00236983, 0.03105614, 0.01542114], ...
                'MoI_at_mass',          [   0.02140958,  0.00027172,  0.00002461; 
                                            0.00027172,  0.01377875, -0.00181920;
                                            0.00002461, -0.00181920,  0.01558906], ...
                'tip_frame',            'F3', ...
                'tip_frame_transformation', struct( ...
                    'q', [ 0 0 0 ], ...
                    'w', [ 1 0 0 ], ...
                    't', pi/2 ...
                ), ...
                "tip_frame_revolute_joint", struct( ...
                    'name', 'J3', ...
                    'axis', [ 0 0 1 ], ... % axis of roation at z-axis
                    'limits', [ -2.8, 2.8 ] ...
                ) ...
            ), ...
            struct( ... 
                'name',                 'arm_B3353', ...
                'base_frame',           'F3', ...
                'mass',                 1.80228141, ...
                'mass_center',          [-0.03825858, 0.20750770, 0.00003309], ...
                'MoI_at_mass',          [    0.05911077, -0.00249612,  0.00000738; 
                                            -0.00249612,  0.00324550, -0.00001767;
                                             0.00000738, -0.00001767,  0.05927043 ], ...
                'tip_frame',            'F4', ...
                'tip_frame_transformation', struct( ...
                    'q', [ common.W_ELBOW_JNT 0 common.L_ARM ], ...
                    'w', [ 1 0 0 ], ...
                    't', -pi/2 ...
                ), ...
                "tip_frame_revolute_joint", struct( ...
                    'name', 'J4', ...
                    'axis', [ 0 0 1 ], ... % axis of roation at z-axis
                    'limits', [ -0.9, 3.1 ] ...
                ) ...
            ), ...
            struct( ... 
                'name',                 'elbow_+_wrist_base_B3355', ...
                'base_frame',           'F4', ...
                'mass',                 2.40016804, ...
                'mass_center',          [0.00498512, -0.00022942, 0.13271662], ...
                'MoI_at_mass',          [    0.01491672,  0.00001741, -0.00150604; 
                                             0.00001741,  0.01482922, -0.00002109;
                                            -0.00150604, -0.00002109,  0.0029446 ], ...
                'tip_frame',            'F5', ...
                'tip_frame_transformation', struct( ...
                    'q', [ -common.W_ELBOW_JNT -common.L_ELBOW 0 ], ...
                    'w', [ 1 0 0 ], ...
                    't', pi/2 ...
                ), ...
                "tip_frame_revolute_joint", struct( ...
                    'name', 'J5', ...
                    'axis', [ 0 0 1 ], ... % axis of roation at z-axis
                    'limits', [ -4.8, 1.3 ] ...
                ) ...
            ), ...
            struct( ... 
                'name',                 'wrist_yaw_B3345', ...
                'base_frame',           'F5', ...
                'mass',                 0.12376019, ...
                'mass_center',          [0.00008921, 0.00511217, 0.00435824], ...
                'MoI_at_mass',          [    0.00005029,  0.00000020, -0.00000005; 
                                             0.00000020,  0.00007582, -0.00000359;
                                            -0.00000005, -0.00000359,  0.0000627 ], ...
                'tip_frame',            'F6', ...
                'tip_frame_transformation', struct( ...
                    'q', [ 0 0 0 ], ...
                    'w', [ 1 0 0 ], ...
                    't', -pi/2 ...
                ), ...
                "tip_frame_revolute_joint", struct( ...
                    'name', 'J6', ...
                    'axis', [ 0 0 1 ], ... % axis of roation at z-axis
                    'limits', [ -1.6, 1.6 ] ...
                ) ...
            ), ...
            struct( ... 
                'name',                 'wrist_pitch', ...
                'base_frame',           'F6', ...
                'mass',                 0.41797364, ...
                'mass_center',          [-0.00012262, -0.01703194, 0.02468336], ...
                'MoI_at_mass',          [    0.00055516, 0.00000061, -0.00000074; 
                                             0.00000061, 0.00024367, -0.00004590;
                                            -0.00000074, -0.00004590, 0.0004535 ], ...
                'tip_frame',            'F7', ...
                'tip_frame_transformation', struct( ...
                    'q', [ 0 0 0 ], ...
                    'w', [ 1 0 0 ], ...
                    't', pi/2 ...
                ), ...
                "tip_frame_revolute_joint", struct( ...
                    'name', 'J7', ...
                    'axis', [ 0 0 1 ], ... % axis of roation at z-axis
                    'limits', [ -2.2, 2.2] ...
                ) ...
            ), ...
            struct( ... 
                'name',                 'wrist_palm', ...
                'base_frame',           'F7', ...
                'mass',                 0.06864753, ...
                'mass_center',          [-0.00007974, 0.00016313, -0.00323552], ...
                'MoI_at_mass',          [    0.00003773, -0.00000019, 0.00000000; 
                                            -0.00000019,  0.00003806, 0.00000000;
                                             0.00000000,  0.00000000, 0.00007408 ], ...
                'tip_frame',            'F8 (Tool Frame)', ...
                'tip_frame_transformation', struct( ... % Tool frame transformation: 
                    'q', [ 0 0 0.06 ], ... % only offset along z-axis
                    'w', [ 0 0 1 ], ... 
                    't', 0 ...
                ), ...
                "tip_frame_revolute_joint", struct( ...
                    'name', 'Unknown', ...
                    'axis', [], ... 
                    'limits', [] ... 
                ) ...
            ) ...
            % struct( ... 
            %     'name',                 'test:wrist_palm+hand+handle', ...
            %     'base_frame',           'F7', ...
            %     'mass',                 0.06864753, ...
            %     'mass_center',          [-0.00007974, 0.00016313, -0.00323552], ...
            %     'MoI_at_mass',          [    0.00003773, -0.00000019, 0.00000000; 
            %                                 -0.00000019,  0.00003806, 0.00000000;
            %                                  0.00000000,  0.00000000, 0.00007408 ], ...
            %     'tip_frame',            'F8 (Tool Frame)', ...
            %     'tip_frame_transformation', struct( ...
            %         'q', [ 0 0 0.06+0.08+0.03 ], ...
            %         'w', [ 0 0 1 ], ...
            %         't', pi ...
            %     ), ...
            %     "tip_frame_revolute_joint", struct( ...
            %         'name', 'J8', ...
            %         'axis', [ 1 0 0 ], ... % axis of roation at z-axis
            %         'limits', [] ... # fixed frame, no joint
            %     ) ...
            % ), ...
            % struct( ... 
            %     'name',                 'test:unknown', ...
            %     'base_frame',           'F8', ...
            %     'mass',                 0.06864753, ...
            %     'mass_center',          [-0.00007974, 0.00016313, -0.00323552], ...
            %     'MoI_at_mass',          [    0.00003773, -0.00000019, 0.00000000; 
            %                                 -0.00000019,  0.00003806, 0.00000000;
            %                                  0.00000000,  0.00000000, 0.00007408 ], ...
            %     'tip_frame',            'F9', ...
            %     'tip_frame_transformation', struct( ...
            %         'q', [ 0 0 0.79 ], ...
            %         'w', [ 1 0 0 ], ...
            %         't', 0 ...
            %     ), ...
            %     "tip_frame_revolute_joint", struct( ...
            %         'name', 'J9', ...
            %         'axis', [ 1 0 0 ], ... % axis of roation at z-axis
            %         'limits', [] ... # fixed frame, no joint
            %     ) ...
            % ), ...
            % struct( ... 
            %     'name',                 'test:unknown', ...
            %     'base_frame',           'F9', ...
            %     'mass',                 0.06864753, ...
            %     'mass_center',          [-0.00007974, 0.00016313, -0.00323552], ...
            %     'MoI_at_mass',          [    0.00003773, -0.00000019, 0.00000000; 
            %                                 -0.00000019,  0.00003806, 0.00000000;
            %                                  0.00000000,  0.00000000, 0.00007408 ], ...
            %     'tip_frame',            'F10', ...
            %     'tip_frame_transformation', struct( ...
            %         'q', [ 0 0 0.84 ], ...
            %         'w', [ 0 0 1 ], ...
            %         't', 0 ...
            %     ), ...
            %     "tip_frame_revolute_joint", struct( ...
            %         'name', 'J10', ...
            %         'axis', [ 0 0 1 ], ... % axis of roation at z-axis
            %         'limits', [] ... # fixed frame, no joint
            %     ) ...
            % ) ...
        ];
    end
    methods(Static)
        % TBD
    end
end
