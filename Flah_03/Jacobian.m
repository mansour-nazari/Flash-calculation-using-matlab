function jacobian = Jacobian(P, T, data, binary, model, nv, sai)
    % Extract required parameters
    z = data(:,1)' / 100;  % Molar composition
    Pc = data(:,2)';        % Critical pressure
    inc = 0.01;            % Increment value
    nc = numel(Pc);        % Number of components
    
    % Preallocate matrices
    jacobian = zeros(nc, nc); 
    K = cell(nc, 2); 
    Phil = cell(nc, 2); 
    Phiv = cell(nc, 2);
    
    % Compute equilibrium ratios and fugacity coefficients
    for i = 1:nc
        perturbation = zeros(1, nc);
        perturbation(i) = 1;  % Perturbation applied to one component at a time
        
        for j = 1:2
            factor = (-1)^(j+1) * inc;
            K{i, j} = exp(sai + factor * perturbation .* sai);
            
            % Compute phase compositions
            liquid_comp = z ./ (1 + nv .* (K{i, j} - 1));
            vapor_comp = K{i, j} .* liquid_comp;
            
            % Compute fugacity coefficients
            [Phi, ~, ~] = EOS(P, T, data, binary, model, liquid_comp, vapor_comp);
            Phil{i, j} = Phi{1};
            Phiv{i, j} = Phi{2};
        end
    end
    
    % Compute Jacobian matrix
    for i = 1:nc
        for j = 1:nc
            numerator = (log(K{j, 1}(i)) - log(Phil{j, 1}(i)) + log(Phiv{j, 1}(i))) - ...
                        (log(K{j, 2}(i)) - log(Phil{j, 2}(i)) + log(Phiv{j, 2}(i)));
            jacobian(i, j) = numerator / (2 * inc * sai(j));
        end
    end
end
