clc
clear
close all
%% compositions
[data,txt] = xlsread('pvt.xlsx','PROP'); %reading components properties
[binary,~] = xlsread('pvt.xlsx','BI');  %reading components BI
CompList = txt(2:end,1);
%% pressure and temp
P = 2500; % Pressure in psia
T = 440+460; % Temperature in R
z = data(:,1)'/100; % mole fraction
%% EOS
EOSList = {'SRK'; 'SRK G&D'; 'PR76'; 'PR78'}; 
[s,~] = listdlg('PromptString','Select an EOS:','SelectionMode','single','ListString', EOSList);
model = EOSList{s};

%% FLASH CALCULATION 
Tc=data(:,3); % CRITICAL TEMP
Pc=data(:,2); % CRITICAL PRESSURE
Omega=data(:,4); % ACENTRIC FACTOR
[K,~, ~] = PhaseStabilityAnalysis(P, T, data, binary, model);
Ki = K; % Initial Equilibrium Ratio Values
sai = log(K); % Creating Sai Numbers
tol = 1; % Total Residual for General While
ct = 0; % General Counter
nv_old = 0.5; %First Guess
alpha = 0.7; % First alpha
error=tol;
IN=1; %indicator
while tol > 1e-7 % General Condition for Sai List
    IN=IN+1;
    if alpha < 0.005 || nv_old > 1 || nv_old < 0
        nv_old = linspace(0.01,0.99,10); % Creating a List for intial Guess
        nv_old = nv_old(randi(10)); % Choosing the initial Guess of the Molar Percentage in Vapour Phase Randomly
        sai = log(Ki); % Assigning the Initial Equilibrium Ratio Values to the  Sai List
        ct = 0; % resetting the Counter
        alpha = 0.7; % resetting the alpha
    else
        sai_old = sai;
        x = z./(1+nv_old*(exp(sai_old)-1)); % The Liqiuid Phase Composition Calculation 
        y = exp(sai_old).*z./(1+nv_old*(exp(sai_old)-1)); % The Vapor Phase Composition Calculation 
        [Phi, Z, Rho] = EOS(P,T,data,binary,model,x,y);
        g = sai_old - log(Phi{1}) + log(Phi{2}); % Creating g Function
        sai = sai_old' - Jacobian(P, T, data,binary, model, nv_old, sai_old)\g'; % Calculation of new Sai Numbers using Jacobian Function-Call
        sai = sai'; % Transposing in order to a row Matrix
        tol = sum(abs(sai - sai_old)); % Assigning the General Residual
        K = exp(sai); % Extracting Equilibrium Ratios from Sai Numbers
        alpha = alpha*0.7;
        ct = ct + 1; % Counter
    end
    nv=rashford_rice(K,z',nv_old);
    nv_old = nv;
    error(IN)=tol;
end

%% IFT
MW=data(:,5); %Molcular weight
pchi=data(:,end); %pachlor number
IFT=surface_tension(x,y,MW',Rho(1),Rho(2),pchi');
%% RESULTS
winprop;
disp(' ')
disp(['EOS: ', model])
disp(['Z: ', num2str(Z(1)), '  , ', num2str(Z(2))])
disp(['Density: ', num2str(Rho(1)), ' lb/cuft  , ', num2str(Rho(2)),' lb/cuft'])
disp(['Nv: ', num2str(nv)])
disp(['IFT: ', num2str(IFT)])

% error table
xx=1:1:length(error);
D.iteration=xx';
D.error=error';

% table of results
struct2table(D)
S.Components = CompList;
S.mole_percent = 100*z';
S.x = 100*x';S.y = 100*y';
% winprop results
if s==1
    S.winporpx=srk(:,1);
    S.winporpy=srk(:,2);
elseif s==2
    S.winporpx=srksg(:,1);
    S.winporpy=srksg(:,2);
    
elseif s==3
    S.winporpx=pr76(:,1);
    S.winporpy=pr76(:,2);
    
elseif s==4
    S.winporpx=pr78(:,1);
    S.winporpy=pr78(:,2);
end
S.PhiL = Phi{1}';
S.PhiV = Phi{2}';S.K = K';
struct2table(S)