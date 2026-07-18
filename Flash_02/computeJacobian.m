function J = computeJacobian(dphi, comp_trial)
    % Builds Jacobian matrix for Newton-Raphson update

    n = length(comp_trial);
    J = zeros(n, n);

    for i = 1:n
        for j = 1:n
            J(i,j) = dphi{2,i}(j) * comp_trial(j);
            if i == j
                J(i,j) = J(i,j) + 1;
            end
        end
    end
end
