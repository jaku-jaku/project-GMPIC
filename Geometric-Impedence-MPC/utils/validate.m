classdef validate
    properties(Constant)
        TOLERANCE = 1e-3;
        VALIDATE_LEVEL_NONE = -1;
        VALIDATE_LEVEL_HARD = 0;
        VALIDATE_LEVEL_SOFT = 1;
    end
    methods(Static)
        function [log_msg, err_msg] = genLogMsg(var_name, true_statement, error_msg)
            log_msg = sprintf("[Lie] [PASS] > validating %s to be %s!", var_name, true_statement);
            err_msg = sprintf("[Lie] [FAIL] > validating %s to be %s, but %s", var_name, true_statement, error_msg);
        end
        function set_validatorLevel(level)
            global gValidatorLevel;
            switch level
                case "all"
                    gValidatorLevel = validate.VALIDATE_LEVEL_SOFT;
                case "basic"
                    gValidatorLevel = validate.VALIDATE_LEVEL_HARD;
                otherwise
                    gValidatorLevel = validate.VALIDATE_LEVEL_NONE;
            end
            helper.loginfo(sprintf("[Lie] > validator level set to %d", gValidatorLevel));
        end
        function judge = compare_float_equal(A, B)
            judge = norm(A-B) < validate.TOLERANCE;
        end
        function if_SO3(var_name, mat_R3x3)
            global gValidatorLevel;
            if gValidatorLevel < validate.VALIDATE_LEVEL_HARD
                return;
            end
            validate.if_dimension(var_name, mat_R3x3, [3,3]);
            validate.if_SO(var_name, mat_R3x3);
            
            dist_to_SO3 = norm(Lie.distance_to_SO3(mat_R3x3));

            % check with verbose:
            [log_msg, err_msg] = validate.genLogMsg(var_name, "SO3", ...
                sprintf("dist_to_SO3=%f",dist_to_SO3));
            assert(dist_to_SO3 < validate.TOLERANCE, err_msg);
            helper.logdebug(log_msg);
        end
        function if_SE3(var_name, mat_R4x4)
            global gValidatorLevel;
            if gValidatorLevel < validate.VALIDATE_LEVEL_HARD
                return;
            end
            mat_R3x3 = mat_R4x4(1:3,1:3);
            
            validate.if_dimension(var_name, mat_R4x4, [4,4]);
            validate.if_SO(var_name, mat_R3x3);
            
            dist_to_SE3 = norm(Lie.distance_to_SE3(mat_R4x4));

            % check with verbose:
            [log_msg, err_msg] = validate.genLogMsg(var_name, "SE3", ...
                sprintf("dist_to_SE3=%f",dist_to_SE3));
            assert(dist_to_SE3 < validate.TOLERANCE, err_msg);
            helper.logdebug(log_msg);
        end
        function if_SO(var_name, mat_R)
            global gValidatorLevel;
            if gValidatorLevel < validate.VALIDATE_LEVEL_HARD
                return;
            end
            if isnumeric(mat_R) 
                % |R| = + 1 %(Murray 2.2)
                det_R = det(mat_R);
                n = length(mat_R);

                % check with verbose:
                [log_msg, err_msg] = validate.genLogMsg(var_name, ...
                    sprintf("in SO{%d}",n), ...
                    sprintf("|R|=%f != 1", det_R));
                assert(validate.compare_float_equal(det_R,1), err_msg);
                helper.logdebug(log_msg);
            end
        end
        function if_so_(var_name, mat_R)
            global gValidatorLevel;
            if gValidatorLevel < validate.VALIDATE_LEVEL_HARD
                return;
            end
            if isnumeric(mat_R) 
                % R.' = - R ===> R.' + R = 0
                delta = mat_R.' + mat_R; %(Murray 2.10)
                n = length(mat_R);
                % check with verbose:
                [log_msg, err_msg] = validate.genLogMsg(var_name, ...
                    sprintf("so{%d}: {R in R^{%dx%d} | R.T + R = 0 }", n), ...
                    sprintf("R.T + R = %f != 0", delta));
                assert(validate.compare_float_equal(det_R,1), err_msg);
                helper.logdebug(log_msg);
            end
        end
        function if_unitVector(var_name, vector)
            global gValidatorLevel;
            if gValidatorLevel < validate.VALIDATE_LEVEL_SOFT
                return;
            end
            % check with verbose:
            [log_msg, err_msg] = validate.genLogMsg(var_name, ...
                "a unit vector", ...
                sprintf("[ %s ]", join(string(vector)) ));
            assert(validate.compare_float_equal(norm(vector), 1), err_msg);
            helper.logdebug(log_msg);
        end
        function if_unitAxis(var_name, w_R3)
            global gValidatorLevel;
            if gValidatorLevel < validate.VALIDATE_LEVEL_SOFT
                return;
            end
            % check with verbose:
            [log_msg, err_msg] = validate.genLogMsg(var_name, ...
                "a unit axis", ...
                sprintf("%s", helper.a2str(var_name,w_R3) ));
            assert(validate.compare_float_equal(det(w_R3), 1), err_msg);
            helper.logdebug(log_msg);
        end
        function if_inRangeStrict(var_name, vector, range)
            global gValidatorLevel;
            if gValidatorLevel < validate.VALIDATE_LEVEL_HARD
                return;
            end
            % check with verbose:
            [log_msg, err_msg] = validate.genLogMsg(var_name, ...
                sprintf(" \in (%f, %f)", range(1), range(2)), ...
                sprintf("%s", helper.a2str(var_name,vector) ));
            assert( range(1) < t_R && t_R < range(2) , err_msg);
            helper.logdebug(log_msg);
        end
        function if_dimension(var_name, mat, dimension)
            global gValidatorLevel;
            if gValidatorLevel < validate.VALIDATE_LEVEL_HARD
                return;
            end
            [n,m] = size(mat);
            % check with verbose:
            [log_msg, err_msg] = validate.genLogMsg(var_name, ...
                sprintf("%dx%d matrix", dimension(1), dimension(2)), ...
                sprintf("it's [%d x %d]", n,m ));
            assert( prod([n,m] == dimension) , err_msg);
            helper.logdebug(log_msg);
        end
        function if_equivalent(A, B, var_name)
            global gValidatorLevel;
            if gValidatorLevel < validate.VALIDATE_LEVEL_HARD
                return;
            end
            % check with verbose:
            [log_msg, err_msg] = validate.genLogMsg(var_name, ...
                sprintf("equivalent!"), ...
                sprintf("not equal :( \n[ Difference ]:\n %s \n", ...
                helper.a2str("A-B",A-B)));
            assert( validate.compare_float_equal(A,B) , err_msg);
            helper.loginfo(log_msg);
        end
    end
end