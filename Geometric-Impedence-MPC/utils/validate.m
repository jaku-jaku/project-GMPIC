classdef validate
    properties(Constant)
        TOLERANCE = 1e-3;
        VALIDATE_LEVEL_NONE = -1;
        VALIDATE_LEVEL_HARD = 0;
        VALIDATE_LEVEL_SOFT = 1;
    end
    methods(Static)
        function message = log(var_name, action)
            message = sprintf("[Lie] > validating %s to be %s", var_name, action);
            helper.logdebug(message);
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
            msg = validate.log(var_name, sprintf("should be SO3! dist_to_SO3=%f",dist_to_SO3));
            assert(dist_to_SO3 < validate.TOLERANCE, msg)
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
            msg = validate.log(var_name, sprintf("should be SE3! dist_to_SE3=%f",dist_to_SE3));
            assert(dist_to_SE3 < validate.TOLERANCE, msg)
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
                msg = validate.log(var_name, ...
                    sprintf("in SO{%d}, with |R|=%f != 1", n, det_R));
                assert(validate.compare_float_equal(det_R,1), msg);
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
                msg = validate.log(var_name, ...
                    sprintf("so{%d}: {R in R^{%dx%d} | R.T + R = 0 } : %f", n, delta));
                assert(validate.compare_float_equal(delta, 0), msg);
            end
        end
        function if_unitVector(var_name, vector)
            global gValidatorLevel;
            if gValidatorLevel < validate.VALIDATE_LEVEL_SOFT
                return;
            end
            msg = validate.log(var_name, sprintf("a unit vector, but: [%s]", join(string(vector))));
            assert(validate.compare_float_equal(norm(vector), 1), msg );
        end
        function if_unitAxis(var_name, w_R3)
            global gValidatorLevel;
            if gValidatorLevel < validate.VALIDATE_LEVEL_SOFT
                return;
            end
            msg = validate.log(var_name, "a unit axis");
            assert(validate.compare_float_equal(det(w_R3), 1), msg );
        end
        function if_inRangeStrict(var_name, vector, range)
            global gValidatorLevel;
            if gValidatorLevel < validate.VALIDATE_LEVEL_HARD
                return;
            end
            msg = sprintf(" \in (%f, %f)", range(1), range(2));
            msg = validate.log(var_name, msg);
            assert( range(1) < t_R && t_R < range(2), msg );
        end
        function if_dimension(var_name, mat, dimension)
            global gValidatorLevel;
            if gValidatorLevel < validate.VALIDATE_LEVEL_HARD
                return;
            end
            [n,m] = size(mat);
            msg = validate.log(var_name, sprintf("%dx%d matrix", dimension(1), dimension(2)));
            assert( prod([n,m] == dimension), msg );
        end
        function test_compare(A, B, name)
            delta = norm(A-B);
            msg = validate.log(name, sprintf("equivalent! epsilon=%f",delta));
            assert( delta < validate.TOLERANCE, msg )
            helper.loginfo("[ x-- test passed! ]");
        end
    end
end