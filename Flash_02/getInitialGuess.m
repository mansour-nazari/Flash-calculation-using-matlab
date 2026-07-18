function logK = getInitialGuess(method, P, T, Tc, Pc, acentric, Mw, mol_frac, data, binary, model)
    % Returns initial logK guess using different methods

    n_comp = length(Pc);

    if method == 1
        % Method 1: Wilson equation
        K = wilsons(T, P, Tc, Pc, acentric);
        logK = log(K);

    elseif method == 2
        % Method 2: Inverse Wilson
        K = wilsons(T, P, Tc, Pc, acentric);
        logK = log(1 ./ K);

    elseif method == 3
        % Method 3: Average of Wilson and its inverse
        K = wilsons(T, P, Tc, Pc, acentric);
        logK = log((K + 1 ./ K) / 2);

    elseif method == 4
        % Method 4: Pure heavy component
        u = ones(1, n_comp) * 1e-5;
        u(Mw == max(Mw)) = 1;
        logK = log(u ./ mol_frac);

    elseif method == 5
        % Method 5: Pure light component
        u = ones(1, n_comp) * 1e-5;
        u(Mw == min(Mw)) = 1;
        logK = log(u ./ mol_frac);

    elseif method == 6
        % Method 6: Fugacity coefficients from EOS
        u0 = zeros(1, n_comp);
        [phi_vals, ~, ~, ~, ~] = EOS_Function(P, T, data, binary, model, mol_frac, u0);
        logK = log(phi_vals{1});
    end
end
