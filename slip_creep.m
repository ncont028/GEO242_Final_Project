function [di_mean, M0_mean, Vd_short, Vd_long] = slip_creep(g_in,alpha,beta)
%slip_creep.m - calculates the slip and creep rates (short and long-term)
%for repeater sequences
%--------------------------------------------------------------------------
% SHORT-TERM CREEP RATE FOR A SEQUENCE
% Chen 2008, eqn. 4
% V_d(i) = [Sum(d(i)) from i=2 to N]/[Sum(Tr(i) frm j=1 to N-1]
% i : individual events in a sequence (family)
% j : recurrence interval in a sequence (family)
% N : number of events in a sequence (family)
% i.e V_d(i) is the short-term creep rate for a sequence (family)
% LONG-TERM CREEP RATE FOR A SEQUENCE
% Chen 2008, eqn. 5
% V_d = [Sum(d(i)) from i=1 to N]/T_obs <-- NOTE THIS TAKES D(I) FOR VERY 1ST
% EVENT!!!!!!!
% T_obs: observation period
%--------------------------------------------------------------------------
% LOAD CORRESPONDING GRID WORKSPACE
g_workspace = ['g_',num2str(g_in),'vars.mat'];
load(g_workspace)

%--------------------------------------------------------------------------
% DEFINE EMPTY CELL ARRAYS
Vd_short = zeros(length(fam_Tr),1); % array to store family (short-term) creep rates
Vd_long = zeros(length(fam_Tr),1); % array to store family (long-term) creep rates
di_mean = zeros(length(fam_Tr),1) ;% slip, di for N events
M0_mean = zeros(length(fam_Tr),1) ;% M0, di for N events

%-------------------------------------------------------------------------
%% FAULT SLIP (di)
%--------------------------------------------------------------------------
% Calculating slip
fam_d = cell(length(fam_Tr),1);
% Loop through families
for f = 1:length(fam_Tr)
    % Loop through events in current family
    for i = 1:length(fam_Tr{f}(:,1))
        M0 = fam_Tr{f}{i,2}(1); % M0 for event i (starts at 2nd event in RE sequence, see Chen 2008)
        d = (10^alpha)*(M0^beta); % slip [cm] for event i (")
        
        % append d and M0 to slip cell array
        fam_d{f,1}{i,1} = d; % cm
        fam_d{f,1}{i,2} = M0; % dyne-cm
    end 
end

%--------------------------------------------------------------------------
%% CALCULATE SHORT- & LONG-TERM CREEP RATES
%--------------------------------------------------------------------------
for f = 1:length(fam_Tr)
    % SHORT-TERM
    d_sum = sum([fam_d{f}{:,1}]); % sum slip for family

    Tr_sum = sum([fam_Tr{f}{:,1}]); % sum Tr for family [sec]
%     Tr_sum = sum([fam_Tr{f}{:,1}]); % sum Tr for family
    Vd_short(f) = d_sum/Tr_sum; % short-term creep rate

    % LONG-TERM
    d1 = (10^alpha)*((fam_Tr{f}{1,2}(2))^beta); % slip for event 1 in a sequence
    d_sum2 = sum([d1,d_sum]); % sum of slip for all events, N, in current family [cm]
    Vd_long(f) = d_sum2/T_obs; % T_obs should be in sec

    % MEAN SLIP & M0 FOR FAMILY
    m0_1 = fam_Tr{f,1}{1,2}(2); % seismic moment for first event
    di_mean(f) = mean([d1,fam_d{f}{:,1}]); % mean slip [cm] for a family (N events)
    M0_mean(f) = mean([m0_1,fam_d{f,1}{:,2}]); % mean m0 [dyne-cm] for a family (N events)

end

%--------------------------------------------------------------------------
%% AVG. CREEP RATE VALUES FOR GRID
%--------------------------------------------------------------------------
% % Once done calculating, find avg. creep rates for the grid
% Vd_short_avg = mean(Vd_short);
% Vd_long_avg = mean(Vd_long);

end