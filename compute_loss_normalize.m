function loss = compute_loss_normalize(X, y, z, gamma, w, nor)
    [m,n] = size(X);
    X_fake = (z*w.' + gamma*[X ones(m, 1)])*inv((w*w.' + gamma*eye(n+1)));
    y_fake = X_fake * w;
    % reverse process
    y_fake = y_fake*nor; y = y*nor;
    % compute loss
    loss = (norm(y_fake - y)^2)/m;
    fprintf('loss=%f\n',loss);
end
