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
IN=1; %indicator
tol = 1; % Total Residual for General While
nv_old = 0.5; %First Guess
error=tol;
while tol > 1e-12 
        IN=IN+1;
        x = z./(1+nv_old*(K-1)); 
        y = K.*z./(1+nv_old*(K-1)); 
        [Phi, Z, Rho] = EOS(P,T,data,binary,model,x,y); 
        K = Phi{1}./Phi{2};
        nv=rashford_rice(K,z',nv_old);
        nv_old = nv; 
        tol = sum((1-((x.*Phi{1})./(y.*Phi{2}))).^2); 
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