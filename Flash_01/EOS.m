function [Phi, Z, Rho] = EOS(P, T, data, binary, model, x, y)
    M = data(:,5)'; % Molecular Weight
    SESRK = data(:,7)'; %volume shift srk
    SEPR = data(:,8)'; % %volume shift PR
    np = 2; % Number of Phases
    Pc = data(:,2)'; % Critical Pressure
    Tc = data(:,3)'+459.67; % Critical Temperature
    af = data(:,4)'; % Acentric Factor
    Phi = cell(1,np); % Pre-allocation for Fugacity Ratios
    Comp = {x,y}; % Pre-allocation for Liquid and Vapor Molar Percent
    Z = cell(1,np); % Pre-allocation for Compressibility Factors
    for i = 1:np
        X = Comp{i};
        R = 10.73159; % Gas Universal Constant
        Tr = T./Tc; % Reduced Temperature Calculation
        switch model
            case 'SRK'  % Soave-Redlich_Kwong
                m = (0.480 + 1.574*af - 0.176*(af.^2));
                aTr = (1 + m.*(1 - sqrt(Tr))).^2;
                delta1 = 0;delta2 = 1;
                omegaa = 0.42748;omegab = 0.08664;
                alpha = (omegaa*(R^2).*(Tc.^2)./Pc);
                a = aTr .* alpha;
                s = sqrt(a).*((X.*sqrt(a))*(1-binary)');
                b = omegab*R*Tc./Pc;
                c = SESRK.*b;
            case 'SRK G&D' % Soave-Redlich_Kwong G&D
                m = (0.48508 + 1.55171*af - 0.15613*(af.^2));
                aTr = (1 + m.*(1 - sqrt(Tr))).^2;
                delta1 = 0;delta2 = 1;
                omegaa = 0.42748;omegab = 0.08664;
                alpha = (omegaa*(R^2).*(Tc.^2)./Pc);
                a = aTr .* alpha;
                s = sqrt(a).*((X.*sqrt(a))*(1-binary)');
                b = omegab*R*Tc./Pc;
                c = SESRK.*b;
            case 'PR76' % Peng-Robinson 76
                m = (0.37464 + 1.54226*af - 0.26992*(af.^2));
                aTr = (1 + m.*(1 - sqrt(Tr))).^2;
                delta1 = 1 - sqrt(2);delta2 = 1 + sqrt(2);
                omegaa = 0.45724;omegab = 0.07780;
                alpha = (omegaa*(R^2).*(Tc.^2)./Pc);
                a = aTr .* alpha;
                s = sqrt(a).*((X.*sqrt(a))*(1-binary)');
                b = omegab*R*Tc./Pc;
                c = SEPR.*b;
            case 'PR78' % Peng-Robinson 78
                m = (0.379642+ 1.48503*af - 0.164423*(af.^2)+ 0.016666*(af.^3));
                aTr = (1 + m.*(1 - sqrt(Tr))).^2;
                delta1 = 1 - sqrt(2); delta2 = 1 + sqrt(2);
                omegaa = 0.457235;omegab = 0.077796;
                alpha = (omegaa*(R^2).*(Tc.^2)./Pc);
                a = aTr .* alpha;
                s = sqrt(a).*((X.*sqrt(a))*(1-binary)');
                b = omegab*R*Tc./Pc;
                c = SEPR.*b;
            otherwise
                error('Unknown EOS, EOS set: [SRK, SRK G&D, PR76, PR78]')    
        end
        at = X*s';
        bt = X*b';
        A = at*P/(R*T)^2;
        B = (bt)*P/(R*T);
        C = -delta2*delta1;
        %% THE COMPRESSIBILTY FACTORS
        p = [1, C*B-1, A-B*(1+C)-(B^2)*(1+2*C), -A*B+C*(B^3+B^2)]; % The Cubic EOS Coefficients
        z = roots(p); % Z Generic EOS Calculation
        if isreal(z) == 1 % Three Real Roots
            z = sort(z);
            zl = z(1); % The Liquid Phase Compressibility Factor
            zv = z(end); % The Vapor Phase Compressibility Factor
        else
            for j = 1:numel(z) % One Real and Two Conjugate Complex Roots
                if isreal(z(j)) == 1
                    zl = z(j);
                    zv = z(j);
                end
            end
        end
        Z{i} = [zl,zv];
        z = Z{i}(i);
        %% FUGACITY COEFFICIENTS
        Phi{i} = exp((b/bt)*(z-1) - log(z-B) - (A/(B*(delta2-delta1)))*(2*s/at-b/bt)*(log((z+delta2*B)/(z+delta1*B)))); % Calculation of the Fugacity Coefficient
    end
    Phi = {Phi{1}.*exp(-P*c/(R*T)), Phi{2}.*exp(-P*c/(R*T))}; % The Liquid Phase and Vapor Fugacity Coefficient
    Z = [Z{1}(1),Z{2}(2)]; % The Liquid and Vapor Compressibility Factor
    Z = [Z(1)*(Z(1)*R*T/P - c*x')/(Z(1)*R*T/P),Z(2)*(Z(2)*R*T/P - c*y')/(Z(2)*R*T/P)]; % The Liquid and Vapor Compressibility Factor Corrected for Volume Shift
    Rho = [(x*M')/(Z(1)*R*T/P), (y*M')/(Z(2)*R*T/P)]; % Liquid and Gas Phase Density lb/cuft
end