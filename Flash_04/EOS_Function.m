function [Phi, dPhidx, dPhidT, Z, Rho] = EOS_Function(P, T, data, binary, model, x, y)
    Pc = data(:,2)'; % Critical Pressure
    Tc = data(:,3)'+459.67; % Critical Temperature
    af = data(:,4)'; % Acentric Factor
    M = data(:,5)'; % Molecular Weight
    SESRK = data(:,7)';
    SEPR = data(:,8)';
    np = 2; % Number of Phases
    nc = numel(Pc); % Number of Components
    Phi = cell(1,np); % Pre-allocation for Fugacity Ratios
    Comp = {x,y}; % Pre-allocation for Liquid and Vapor Molar Percent
    Z = cell(1,np); % Pre-allocation for Compressibility Factors
    dPhidx = cell(np,nc); % Pre-allocation for Fugacity Derivatives to X
    dPhidT = cell(np,1); % Pre-allocation for Fugacity Derivatives to T
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
        u = C+1;
        w = -C;
        dradadT = -0.5*sqrt(alpha).*m.*(1./Tc).*(Tr.^(-0.5));
        dsdT = dradadT.*((X.*sqrt(a))*(1-binary)') + sqrt(a).*((X.*dradadT)*(1-binary)');
        datdT = X*dsdT';
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
        %% FUGACITY COEFFICIENTS AND DERIVAITVES
        Phi{i} = exp((b/bt)*(z-1) - log(z-B) - (A/(B*(delta2-delta1)))*(2*s/at-b/bt)*(log((z+delta2*B)/(z+delta1*B)))); % Calculation of the Fugacity Coefficient
        for j = 1:nc % Creating Loop for the Calculation of the Fugacity Derivatives
            for k = 1:nc
                dbdx = b(k);
                dadx = 2*s(k);
                dBdx = B*dbdx/bt;
                dAdx = A*dadx/at;
                r = -A*dBdx - B*dAdx + 2*C*B*dBdx + 3*C*(B^2)*dBdx;
                dzdx = -((C*dBdx)*(z^2) + (dAdx-(1+C)*dBdx - 2*B*(1+2*C)*dBdx)*z + r)/(3*(z^2) + 2*(C*B-1)*z + A-(1+C)*B - (1+2*C)*(B^2)) - P*c(k)/(R*T);
                term1 = (b(j)/bt)*dzdx - (b(j)/bt^2)*(z-1)*dbdx;
                term2 = - (dzdx-dBdx)/(z-B);
                term31 = dAdx/B - (A*dBdx)/(B^2);
                term32 = log((2*z + B*(u+sqrt(u^2-4*w)))/(2*z + B*(u-sqrt(u^2-4*w))));
                term33 = (2*dzdx+dBdx*(u+sqrt(u^2-4*w)))/((2*z + B*(u+sqrt(u^2-4*w))));
                term34 = (2*dzdx+dBdx*(u-sqrt(u^2-4*w)))/((2*z + B*(u-sqrt(u^2-4*w))));
                term3 = (b(j)/bt-(2*s(j))/at)*(term31*term32+(A/B)*(term33-term34))/sqrt(u^2-4*w);
                term41 = -b(j)*dbdx/(bt^2) + dadx*(2*s(j))/(at^2) - 2*(1-binary(j,k))*sqrt(a(j)*a(k))/at;
                term4 = A*term32*term41/(B*sqrt(u^2-4*w));
                dPhidx{i,j}(k) = term1+term2+term3+term4;
            end
        end
        dAdT = datdT*A/at - 2*A/T;
        dBdT = -B/T;
        r = -A*dBdT - B*dAdT + 2*C*B*dBdT + 3*C*(B^2)*dBdT;
        dzdT = -((C*dBdT)*(z^2) + (dAdT-(1+C)*dBdT - 2*B*(1+2*C)*dBdT)*z + r)/(3*(z^2) + 2*(C*B-1)*z + A-(1+C)*B - (1+2*C)*(B^2))-(X*c')/(R*T);
        term1 = (b/bt)*dzdT - (dzdT-dBdT)/(z-B);
        term2 = (b/bt - 2*s/at)*(dAdT/B-A*dBdT/B^2)*(log((2*z+B*(u+sqrt(u^2-4*w)))/(2*z+B*(u-sqrt(u^2-4*w)))));
        term3 = (b/bt - 2*s/at)*(A/B)*(((2*dzdT+dBdT*(u+sqrt(u^2-4*w)))/(2*z+B*(u+sqrt(u^2-4*w))))-((2*dzdT+dBdT*(u-sqrt(u^2-4*w)))/(2*z+B*(u-sqrt(u^2-4*w)))));
        term4 = ((2*s/at^2)*(datdT)-(2/at)*(dsdT))*(A/B)*(log((2*z+B*(u+sqrt(u^2-4*w)))/(2*z+B*(u-sqrt(u^2-4*w)))));
        dPhidT{i} = term1 + (1/sqrt(u^2-4*w))*(term2+term3+term4);
    end
    Phi = {Phi{1}.*exp(-P*c/(R*T)), Phi{2}.*exp(-P*c/(R*T))}; % The Liquid Phase and Vapor Fugacity Coefficient
    Z = [Z{1}(1),Z{2}(2)]; % The Liquid and Vapor Compressibility Factor
    Z = [Z(1)*(Z(1)*R*T/P - c*x')/(Z(1)*R*T/P),Z(2)*(Z(2)*R*T/P - c*y')/(Z(2)*R*T/P)]; % The Liquid and Vapor Compressibility Factor Corrected for Volume Shift
    Rho = [(x*M')/(Z(1)*R*T/P), (y*M')/(Z(2)*R*T/P)]; % Liquid and Gas Phase Density lb/cuft
end