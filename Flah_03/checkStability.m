function [num_phases, K_final] = checkStability(logK, comp_trial, mol_frac, Mw, P, T, R, data, binary, model)
    % Determines whether the system is stable or splits into two phases

    % Check trivial stability (composition sum)
    if sum(comp_trial) <= 1 + 1e5 * eps
        num_phases = 1;
        K_final = zeros(1, length(mol_frac));
        return
    end

    % Get Z-factors for feed and trial phases
    [~, ~, ~, Zf, ~] = EOS_Function(P, T, data, binary, model, mol_frac, comp_trial);

    % Densities for comparison
    rho_feed  = P * dot(mol_frac, Mw) / (Zf(1) * R * T);
    rho_trial = P * dot(comp_trial, Mw) / (Zf(2) * R * T);

    % Choose proper K-values based on phase stability
    if rho_feed > rho_trial
        K_final = exp(logK);
    else
        K_final = 1 ./ exp(logK);
    end

    num_phases = 2;
end
