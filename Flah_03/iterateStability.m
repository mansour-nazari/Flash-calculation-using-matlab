function [logK, comp_trial, converged] = iterateStability(logK_init, mol_frac, data, binary, model, P, T)
    % Performs Newton-Raphson iteration to find trial composition
    % Returns converged logK, trial composition, and status

    logK = logK_init;
    comp_trial = mol_frac .* exp(logK);
    max_iter = 200;
    tolerance = 1e-5;
    iter = 0;
    converged = true;

    while iter < max_iter
        comp_trial = comp_trial / sum(comp_trial);
        logK_old = logK;

        % Get fugacity coefficients and derivatives
        [phi_vals, dphi, ~, ~, ~] = EOS_Function(P, T, data, binary, model, mol_frac, comp_trial);

        % Residual vector for NR method
        g = logK_old - log(phi_vals{1}) + log(phi_vals{2});

        % Build Jacobian matrix
        J = computeJacobian(dphi, comp_trial);

        % Newton update
        delta = J \ g';
        logK = logK_old' - delta;
        logK = logK';

        % Update trial composition
        comp_trial = mol_frac .* exp(logK);

        % Convergence check
        if sum(abs(logK - logK_old)) < tolerance
            return
        end
        iter = iter + 1;
    end

    % Did not converge within limit
    converged = false;
end
