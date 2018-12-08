function [uOptReshape, rMax] = OptimizeNonlinear()

    %OptimizeLinear Performs nonlinear optimization of the constellation
    %separation problem
    %   Performs nonlinear optimization of the constellation separation problem

    %% Import global variables

    global Amin;
    global Amax;
    global N;
    global dt;
    global D;
    global delta_des;
    global epsTheta;
    global epsOmega;
    global re;
    global T;
    global r0;
    global w0;
    global theta0;

    %% Create cost function
    f = [zeros(1,N*T), 1];
    x0 = [repmat(Amin+0.01, N*T, 1); -(475e3+re)];  % Initial guess

    costfun = @(x)f*x;
    
    %Create (highly) nonlinear constraint function

    function [c, ceq] = nonlcon(x)
  
        t = x(end);
        
        %Reshape x into form suitable for computing trajectory
        u = reshape(x(1:end-1), T, N);
        u = u';
        
        %Calculate state at time T
        [r, w, theta] = trajectory(u);
        rT = r(:,end);
        wT = w(:,end);
        thetaT = theta(:,end);
        
        %Compute constraints on r, theta, and w
        rConst = max(-rT) - t;
        thetaConst = max(D*thetaT - delta_des) - epsTheta;
        wConst = max(D*wT) - epsOmega;
        
        c = [rConst; thetaConst; wConst];
        ceq = [];
        
    end

    %Formulate other constraints, mostly empty
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [repmat(Amin, N*T, 1); -Inf];
    ub = [repmat(Amax, N*T, 1); Inf];
    
    
    % Perform optimization. NOTE: for final results, disable parallelization
    options = optimoptions('fmincon', 'Display', 'iter', 'MaxFunctionEvaluations', 1000*N*T, 'UseParallel', false);
    result = fmincon(costfun, x0, A, b, Aeq, beq, lb, ub, @nonlcon, options);

    uOpt = result(1:end-1);  % Extract area commands
    thresh = result(end);  % Extract radius

    uOptReshape = reshape(uOpt, T, N);
    uOptReshape = uOptReshape';

    rMax = -thresh;  % Determines maximized radius

end















