classdef validate
    methods(Static)
        function message = log(var_name, action)
            message = sprintf("[Lie] > validating %s to be %s", var_name, action);
            helper.logdebug(message);
        end
        function set_validatorLevel(level)
            global gValidatorLevel;
            switch level
                case "all"
                    gValidatorLevel = 1;
                case "basic"
                    gValidatorLevel = 0;
                otherwise
                    gValidatorLevel = -1;
            end
            helper.loginfo(sprintf("[Lie] > validator level set to %d", gValidatorLevel));
        end
        function if_SO3(var_name, mat_R3x3)
            global gValidatorLevel;
            if gValidatorLevel < 0
                return;
            end
            validate.if_dimension(var_name, mat_R3x3, [3,3]);
            validate.if_SO(var_name, mat_R3x3);
        end
        function if_SE3(var_name, mat_R4x4)
            global gValidatorLevel;
            if gValidatorLevel < 0
                return;
            end
            mat_R3x3 = mat_R4x4(1:3,1:3);
            
            validate.if_dimension(var_name, mat_R4x4, [4,4]);
            validate.if_SO(var_name, mat_R3x3);
            msg = validate.log(var_name, "should be homogeneous!");
            assert(mat_R4x4(4,4) == 1, msg)
        end
        function if_SO(var_name, mat_R)
            global gValidatorLevel;
            if gValidatorLevel < 0
                return;
            end
            if isnumeric(mat_R) 
                % |R| = + 1, RR' = I
                delta = abs(det(mat_R) - 1); %(Murray 2.2)
                n = length(mat_R);
                msg = validate.log(var_name, ...
                    sprintf("in SO{%d}, with |R|=1, |R|-1=%f", n, delta));
                assert(delta < 1e-8, msg);
            end
        end
        function if_so_(var_name, mat_R)
            global gValidatorLevel;
            if gValidatorLevel < 0
                return;
            end
            if isnumeric(mat_R) 
                % R.' = - R ===> R.' + R = 0
                delta = mat_R.' + mat_R; %(Murray 2.10)
                n = length(mat_R);
                msg = validate.log(var_name, ...
                    sprintf("so{%d}: {R in R^{%dx%d} | R.T + R = 0 } : %f", n, delta));
                assert(delta < 1e-8, msg);
            end
        end
        function if_unitVector(var_name, vector)
            global gValidatorLevel;
            if gValidatorLevel < 0
                return;
            end
            msg = validate.log(var_name, "a unit vector");
            assert( norm(vector) == 1, msg );
        end
        function if_unitAxis(var_name, w_R3)
            global gValidatorLevel;
            if gValidatorLevel < 0
                return;
            end
            msg = validate.log(var_name, "a unit axis");
            assert( det(w_R3) == 1, msg );
        end
        function if_inRangeStrict(var_name, vector, range)
            global gValidatorLevel;
            if gValidatorLevel < 0
                return;
            end
            msg = sprintf(" \in (%f, %f)", range(1), range(2));
            msg = validate.log(var_name, msg);
            assert( range(1) < t_R && t_R < range(2), msg );
        end
        function if_dimension(var_name, mat, dimension)
            global gValidatorLevel;
            if gValidatorLevel < 0
                return;
            end
            [n,m] = size(mat);
            msg = validate.log(var_name, sprintf("%dx%d matrix", dimension(1), dimension(2)));
            assert( prod([n,m] == dimension), msg );
        end
    end
end