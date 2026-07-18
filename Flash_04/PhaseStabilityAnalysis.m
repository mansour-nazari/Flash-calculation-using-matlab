function [K_final, num_phases, comp_trial] = PhaseStabilityAnalysis(P, T, data, binary, model)

    % -----------------------------------------------------------
    % PHASE STABILITY ANALYSIS USING TPD METHOD
    % -----------------------------------------------------------
    % Inputs:
    %   P        - Pressure [psia]
    %   T        - Temperature [F]
    %   data     - Component data matrix
    %   binary   - Binary interaction parameters
    %   model    - EOS model type
    % Outputs:
    %   K_final      - Final K values
    %   num_phases   - Number of phases (1 or 2)
    %   comp_trial   - Trial composition
    % -----------------------------------------------------------
    
    
    % Constants & Inputs
    R        = 10.7335;                        % Gas constant (psia.ft3/lbmol.R)
    mol_frac = data(:,1)' / 100;                % Molar composition (converted to fraction)
    Pc       = data(:,2)';                      % Critical pressures
    Tc       = data(:,3)' + 459.67;             % Critical temperatures (R)
    acentric = data(:,4)';                      % Acentric factors
    Mw       = data(:,5)';                      % Molecular weights
    n_comp   = length(Pc);                      % Number of components

    % Try different initial guesses for logK
    for method = 1:6
        % Generate initial guess
        logK = getInitialGuess(method, P, T, Tc, Pc, acentric, Mw, mol_frac, data, binary, model);
        
        % Run iteration
        [logK_final, comp_trial, converged] = iterateStability(logK, mol_frac, data, binary, model, P, T);

        % Skip if not converged
        if ~converged
            continue
        end

        % Check if system is stable or splits into two phases
        [num_phases, K_final] = checkStability(logK_final, comp_trial, mol_frac, Mw, P, T, R, data, binary, model);
        
        if num_phases == 2
            break
        end
    end
end
