function [jacobian] = JacobianCalc(P, T, data, binary, model, nv, sai)
    % Extract data
    z = data(:,1)' / 100;  % Molar composition
    Pc = data(:,2)';       % Critical pressure
    nc = numel(Pc);        % Number of components
    
    % Compute equilibrium ratios and phase compositions
    K = exp(sai);
    x = z ./ (1 + nv * (K - 1));
    y = K .* z ./ (1 + nv * (K - 1));
    
    % Compute fugacity coefficient derivatives
    [~, dPhidx, ~] = EOS_Function(P, T, data, binary, model, x, y);
    
    % Initialize Jacobian matrix
    jacobian = zeros(nc, nc);
    
    % Compute Jacobian elements
    for i = 1:nc
        for j = 1:nc
            term1 = -dPhidx{1, i}(j) * K(j) * ((-nv * z(j)) / (1 + nv * (K(j) - 1))^2);
            term2 = dPhidx{2, i}(j) * y(j) * ((1 - nv) / (1 + nv * (K(j) - 1)));
            jacobian(i, j) = term1 + term2;
            
            % Add identity matrix component
            if i == j
                jacobian(i, j) = jacobian(i, j) + 1;
            end
        end
    end
end