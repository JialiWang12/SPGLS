function w_opt = rank_decompose(M, W)
    n = size(W, 1);
    q_forms = ones(n, 1);
    X_new = W;
    [V, D] = eig(X_new);
    u = zeros(n, n);
    
    for i =1:n
        u(:, i) = V(:, i)*sqrt(D(i, i));
        q_forms(i) = u(:, i).'*M*u(:, i);
    end
    
    while(~(q_forms <= 0))
        pos_idx = -1;
        neg_idx = -1;
        
        for i = 1:n
            if(q_forms(i) > 0 &&  pos_idx == -1)
                pos_idx = i;
            elseif(q_forms(i) < 0 && neg_idx == -1)
                neg_idx = i;
            end      
        end
        
        a = q_forms(pos_idx);
        b = 2*u(:, pos_idx).'*M*u(:, neg_idx);
        c = q_forms(neg_idx);
        t = (-b + sqrt(b*b - 4*a*c))/(2*a);
        u_pos = (t / sqrt(t*t + 1))*u(pos_idx) + (1 / sqrt(t*t + 1))*u(neg_idx);
        u_neg = - (1 / sqrt(t*t + 1))*u(pos_idx) + (t / sqrt(t*t + 1))*u(neg_idx);
        u(pos_idx) = u_pos;
        u(neg_idx) = u_neg;
        q_forms(pos_idx) = u(pos_idx).'*M*u(pos_idx);
        q_forms(neg_idx) = u(neg_idx).'*M*u(neg_idx);
    end  
    w_opt = u(:, end);
end