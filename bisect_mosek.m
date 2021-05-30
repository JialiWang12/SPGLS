function [w_opt, W_opt, B, Fq, total_time, iter] = bisect_mosek(X,y, z, gamma, bias_upper, q_start, tol,varargin)


    a = 0.0;
    b = q_start;
    time_list = [];
    iter = 0;
    [opt, w_opt, W_opt, time] = evaluate_Fq(X, y, z, b, gamma, bias_upper);
    time_list = [time_list time];
    iter = iter + 1;
    Bq = build_quad_B(X, y, z, b, gamma);
    bb = build_lin_B(X, y, z, b, gamma);
    B = [Bq, bb; bb.', 0];
    
    
    while (b -a) > tol
        q_test = (a+b) / 2;
        [Fq, w_q, W, time] = evaluate_Fq(X, y, z, q_test, gamma, bias_upper);
        time_list = [time_list time];
        iter = iter + 1;     
        if Fq >= 0
            a = q_test;
        else
            b = q_test;
            w_opt = w_q;
            W_opt = W;
            Bq = build_quad_B(X, y, z, q_test, gamma);
            bb = build_lin_B(X, y, z, q_test, gamma);
            B = [Bq, bb; bb.', 0];
            
%             fval(X,y,z,gamma,w_opt(1:12))
        end
        
    end
   total_time = sum(time_list);
end

function A = build_A(X, y, z, q, gamma)
    m = size(X, 1);
    n = size(X, 2);
    A_11 = X.'*X;
    A_12 = X.'*ones(m, 1);
    A_13 = (1.0 / gamma)*X.'*(z - y);
    A_14 = -X.'*y;
    
    A_1 = [A_11, A_12, A_13, A_14];
    
    A_21 = ones(m, 1).'*X;
    A_22 = ones(m, 1).'*ones(m, 1);
    A_23 = (1.0 / gamma) * ones(m, 1).'*(z - y);
    A_24 = -ones(m, 1).'*y;
    
    A_2 = [A_21, A_22, A_23, A_24];
    
    A_31 = (1.0 / gamma) * (z - y).'*X;
    A_32 = (1.0 / gamma) * (z-y).'*ones(m, 1);
    A_33 = (1.0 / gamma) * (1.0 / gamma) * (z-y).'*(z-y) - (1.0 / gamma) * (1.0 / gamma) * q;
    A_34 = (1.0 / gamma) * y.'*(y - z) - (1.0 / gamma) * q;
    
    A_3 = [A_31, A_32, A_33, A_34];
    
    A_41 = -y.'*X;
    A_42 = -y.'*ones(m, 1);
    A_43 = (1.0 / gamma) * y.'*(y - z) - (1.0 / gamma)* q;
    A_44 = y.'*y - q;
    
    A_4 = [A_41, A_42, A_43, A_44];
    
    A = [A_1; A_2; A_3; A_4];
end

function B = build_B(X, y, z, q, gamma)
    m = size(X, 1);
    n = size(X, 2);
    
    B_11 = eye(n);
    B_12 = zeros(n, 3);
    
    B_1 = [B_11, B_12];
    
    B_2 = zeros(3, n+3);
    B_2(2, n+3) = -0.5;
    B_2(3, n+2) = -0.5;
    
    B = [B_1; B_2];
end

function C = build_C(X, y, z, q, gamma)
    n = size(X, 2);
    C = zeros(n+3, n+3);
    C(n+3, n+3) = -1.0;
end

function Aq = build_quad_A(X, y , z , q, gamma)
    m = size(X, 1);
    n = size(X, 2);

    Aq_11 = (1.0 / gamma) * (1.0 / gamma) * (z-y).'*(z-y) - (1.0 / gamma) * (1.0 / gamma) * q;
    Aq_12 = (1.0 / gamma) * (z - y).'*X;
    Aq_13 = (1.0 / gamma) * (z-y).'*ones(m, 1);
    
    Aq_1 = [Aq_11, Aq_12, Aq_13];
    
    Aq_21 = (1.0 / gamma)*X.'*(z - y);
    Aq_22 = X.'*X;
    Aq_23 = X.'*ones(m, 1);
    
    Aq_2 = [Aq_21, Aq_22, Aq_23];
    
    Aq_31 = (1.0 / gamma) * ones(m, 1).'*(z - y);
    Aq_32 = ones(m, 1).'*X;
    Aq_33 = ones(m, 1).'*ones(m, 1);
    
    Aq_3 = [Aq_31, Aq_32, Aq_33];
    
    Aq = [Aq_1; Aq_2; Aq_3];
end

function aq = build_lin_A(X, y , z, q, gamma)
    m = size(X, 1);
    n = size(X, 2);
    
    aq_1 = (1.0 / gamma) * y.'*(y - z) - (1.0 / gamma) * q;
    aq_2 = -X.'*y;
    aq_3 = -ones(m, 1).'*y;
    
    aq = [aq_1; aq_2; aq_3];
end

function ac = build_const_A(X, y, z, q, gamma)
    ac = y.'*y - q;
end

function Bq = build_quad_B(X, y, z, q, gamma)
    m = size(X, 1);
    n = size(X, 2);
    
    Bq = eye(n+2);
    Bq(1, 1) = 0;
end

function bb = build_lin_B(X, y, z, q, gamma)
    m = size(X, 1);
    n = size(X, 2);
    
    bb = zeros(n+2, 1);
    bb(1) = -1;
end

function ub = bias_bound_upper(X, y, z, q, gamma)
    m = size(X, 1);
    n = size(X, 2);
    
    ub = zeros(n+2, 1);
    ub(n+2) = 1.0;
end

function lb = bias_bound_lower(X, y, z, q, gamma)
    m = size(X, 1);
    n = size(X, 2);
    
    lb = zeros(n+2, 1);
     lb(n+2) = -1.0;
end

function Mb = bias_matrix_bound(X, y, z, q, gamma)
    m = size(X, 1);
    n = size(X, 2);
    
    Mb = zeros(n+2, n+2);
    Mb(n+2, n+2) = 1.0;
end


function [val_opt, w_opt, W_opt,time] = evaluate_Fq(X, y, z, q, gamma, bias_upper)
    Aq = build_quad_A(X, y, z , q, gamma);
    ab = build_lin_A(X, y, z, q, gamma);
    ac = build_const_A(X, y, z, q, gamma);

    Bq = build_quad_B(X, y, z, q, gamma);
    bb = build_lin_B(X, y, z, q, gamma);

    ub = bias_bound_upper(X, y, z, q, gamma);
    lb = bias_bound_lower(X, y, z, q, gamma);
    Mb = bias_matrix_bound(X, y, z, q, gamma);
    
    m = size(X, 1);
    n = size(X, 2);
    
    % begin mosek
    temp0 = tic;
    clear prob
    [r, res] = mosekopt('symbcon');
    A_bar = [Aq, ab; ab', ac];
    B1_bar = [Bq, 0.5*bb; 0.5*bb',0];
    B2_bar = zeros(n+3,n+3); B2_bar(end,end) = 1;
    % coefficients in the objective
    [r0 c0 v0] = find(tril(A_bar));
    % coefficients in the constrains
    [r1 c1 v1] = find(tril(B1_bar));
    [r2 c2 v2] = find(tril(B2_bar));
    % idx in the objective
    prob.bardim = [n+3,n+3];
    prob.barc.subj = ones(length(v0),1);
    prob.barc.subk = [r0];
    prob.barc.subl = [c0];
    prob.barc.val = [v0];
    % idx in the constrains
    prob.blc = [0,1];
    prob.buc = [0,1];
    prob.a = sparse([], [], [], 2, 0);
    prob.bara.subi = [ones(length(v1),1);2];
    prob.bara.subj = [ones(1+length(v1),1)];
    prob.bara.subk = [r1;r2];
    prob.bara.subl = [c1;c2];
    prob.bara.val = [v1;v2];

    param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1.0e-12;
    [r,res] = mosekopt('minimize info echo(0)',prob, param);
    ttotal = toc(temp0); time = ttotal;
    val_opt = res.info.MSK_DINF_INTPNT_PRIMAL_OBJ;
    W_bar = zeros(n+3);
    [r0 c0 ~ ]= find(tril(ones(n+3)));
    subsrj =  sub2ind([n+3,n+3],r0,c0);W_bar(subsrj) =  res.sol.itr.barx(1:length(subsrj));
    W_bar = W_bar + tril(W_bar,-1)';
    W_opt = W_bar(1:n+2,1:n+2); w_opt = W_bar(1:n+2,end);

end
