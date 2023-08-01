classdef Error
    properties
        chi
    end
    properties (Constant)
        chiRef = reshape(rand(4,4), [4,4]);
        Dm  = eye(6);
        Q =0.2*eye(6);
        R = 10*eye(6);
    end
    methods (Static)
        function [psiX, psiXse] = ConfgError(chi)
            chi = reshape(chi, [4,4]);
            psiX = trace(eye(3) - Error.chiRef(1:3,1:3)'*chi(1:3,1:3)) + 1/2 * (chi(1:3,4) - Error.chiRef(1:3,4))'*(chi(1:3,4) - Error.chiRef(1:3,4));
            psiXse = logm(psiX);
        end

        function [U,err,X] = eMPC(N,xi0,xiRef, dxiRef, mat, Mr, Gr)
             import casadi.*
             % creat casadi obj & variables
            opti = casadi.Opti();
            dpsi = opti.variable(6,N+1);
            u = opti.variable(6,N+1);
            xi = opti.variable(6,N+1);
            cost = 0;
            for r = 1:N-1
                 cost = cost + xi(:,r)'*Error.Q*xi(:,r) + dpsi(:,r)'*Error.Q*dpsi(:,r) + u(:,r)'*Error.R*u(:,r);
            end
            [Qf, ~, ~] = dare(mat, xiRef, Error.Q,1,[],[], 'anti');
            TerminalCost = dpsi(:,N)'*Qf*dpsi(:,N);
            % cost-to-go + terminal cost
            opti.minimize(cost + TerminalCost);

            for k = 1:N-2
                 opti.subject_to(dpsi(:,k+1) == -mat*dpsi(:,k) + xiRef-xi0)
                opti.subject_to(xi(:,k) ==Error. Dm\Mr*dxiRef + xiRef  + Gr - u(:,k))
             end
            opti.subject_to(dpsi(:,N) == zeros(6,1)) % terminal velocity equals zero
 
            % solver options
            opts = struct;
            opts.ipopt.max_iter = 100;
            opts.ipopt.print_level = 3; %0,3
            opts.print_time = 1; %0,1
            opts.ipopt.acceptable_obj_change_tol = 1e-6;
            opti.solver('ipopt', opts);
            sol = opti.solve();
            U = sol.value(u);
            err = sol.value(dpsi);
            X = sol.value(xi);

        end
    end
end
