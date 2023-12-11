%% SET UP CODE ============================================================

global barD  % Bar depending on the operative system to acces folders
format     long 
fs         = 22;        % fontsize value
lw         = 1.2;       % linewith value
figuresDir = 'figures'; % name of the general folder for figures
pwdDir     = pwd;       % actual directory to make the references
warning('off','all')
beep off

% Create folders for plots ------------------------------------------------
mkdir data
mkdir figures
mkdir figures_doubleGyre
mkdir figures_UF
mkdir figures_vKVS
mkdir figures_sine
% mkdir figures_MC_1D
% mkdir figures_MC_2D
% mkdir figures_MC_GCI
% mkdir figures_MCvsMoD
% mkdir figures_MoD
% mkdir figures_MoD_1D
% mkdir figures_paper
% mkdir figures_runMoD

% This avoid problems with saving figures----------------------------------
set(0, 'DefaultFigureRenderer', 'painters');

% This detects if you are in Windows or Linux to change the bar / or \ ----
for i=length(pwdDir):-1:1 % to find the folder name and if windows/linux
    if pwdDir(i)=='/'     % linux system = 1
        barD=pwdDir(i);
        system=1;
        ks=i;
        break
    elseif pwdDir(i)=='\' % windows system = 0
        barD=pwdDir(i);
        system=0;
        ks=i;
        break
    end
end
nameDir = pwdDir(ks+1:end); % Current folder name

% Adding paths ------------------------------------------------------------
addpath([pwd barD 'matlab_functions'])

% cd ..
% cd('solver')
% here=pwd;
% addpath(here) 
% cd ..
% cd(pwdDir)

% Colors to use in plots --------------------------------------------------
color1 = [0.0000 0.4470 0.7410];
color2 = [0.8500 0.3250 0.0980];
color3 = [0.9290 0.6940 0.1250];
color4 = [0.4940 0.1840 0.5560];
color5 = [0.4660 0.6740 0.1880];
color6 = [0.3010 0.7450 0.9330];
color7 = [0.6350 0.0780 0.1840];
color8  = [0    0    1   ]; % intense blue
color9  = [0    0.5  0   ]; % dark green
color10 = [1    0    0   ]; % intense red
color11 = [0    0.75 0.75]; % cian
color12 = [0.75 0    0.75]; 
color13 = [0.75 0.75 0   ];
color14 = [0.25 0.25 0.25]*0; % black
color15 = [0.5 0.5 0.5]*1; % grey



