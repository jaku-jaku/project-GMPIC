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
            if isnumeric(mat_R3x3) 
                delta = abs(det(mat_R3x3) - 1);
                msg = validate.log(var_name, sprintf("in SO3, with |R|=1, |R|-1=%f", delta));
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
        function if_dimension(var_name, mat, dimension)
            global gValidatorLevel;
            if gValidatorLevel < 0
                return;
            end
            msg = validate.log(var_name, sprintf("%dx%d matrix", dimension(1), dimension(2)));
            assert( norm(size(mat) - dimension) == 0, msg );
        end
    end
end