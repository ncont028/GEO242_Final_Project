%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAC_final_project.m
% NAC
% GEO 242 (F22)
% Calculate fault slip & creep for the Hayward fault
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEFORE RUNNING THIS SCRIPT:
% 1) run find_empty_lines.sh
    % Find "empty"* lines in files.
    % https://www.quora.com/How-do-you-list-only-the-empty-lines-in-a-file-using-grep
    % *Finds lines without ":" (which includes header):
    % grep -vn : RE_catalog_001.txt | awk -F: '{print $1}'>>g1_empty.csv"
% 2) Copy RE catalog files & empty csv files into:
    % 'Hayward_RE_catalog" dir containing RE catalogs
    % 'HF_empty' dir containing empty lines csv files
% 3) Copy fault grids.txt into current working directory
%--------------------------------------------------------------------------
% SECTIONS:
% 1) Loading RE catalogs, calculating Tr & M0
% 2) beta_d calibrations + jackknife sampling
% 3) USGS creep rate data
% 4) alpha_estimates
% 5) slip and creep rate calculations for RE grids
% ALL TIMES ARE IN SEC!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET UP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD IN RE CATALOG FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sample:
% start time                lon       lat      depth  mag
% 1994-05-05 10:41:56.7200 -122.474 38.1567 6.8 1.64 401696
g_files = ls('Hayward_RE_catalogs/RE_catalog_*.txt'); % file location
[num_gs blah] = size(g_files);

% # empty lines = number of families
g_empty_files = ls('HF_empty/g*_empty.csv'); % empty lines files; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETTING UP VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FAMILY SPLITING
val_st = ''; % for strings
val_n = nan; % for numbers
% RECURRENCE INTERVALS
Tr_min = (0.1)*(60*60*24*365); % minimum recurrence interval; 0.1 yrs, as defined by Chen et al (2008).
                               % yr --> sec
% SLIP CALCULATIONS
% Parkfield parameters (Nadeau 1998)
alphaP = -2.36;
betaP = 0.17;

% Khoshmanesh (2015): α = −1.56, β = 0.10
alphaK = -1.56; betaK = 0.10;

% CREEP RATES
T_obs = (2020-1984)*(60*60*24*365); % observation period [yrs]

%***********************
% Spring cleaning
%***********************
clear blah
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATING RECURRENCE INTERVALS, SEISMIC MOMENTS, ETC.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOP THROUGH EACH GRID

for g=1:num_gs
    %----------------------------------------------------------------------
    % LOAD RE CATALOG & EMPTY LINES FILE
    %----------------------------------------------------------------------
    % open RE catalog
    fileID = fopen(['Hayward_RE_catalogs/',g_files(g,:)]);
    g_file = textscan(fileID,'%s %s %f %f %f %f %d','HeaderLines',1);
    fclose(fileID);

    % load empty lines csv
    empty_str = g_empty_files(g,:);
    g_empty = load(['HF_empty/',empty_str])-1; % subtract by 1 to not count header

    %----------------------------------------------------------------------
    % SEPARATE CATALOG INTO FAMILIES
    %----------------------------------------------------------------------
    % length of g_empty = number of families
    for n=1:length(g_empty)
        % index of new value
        idx = g_empty(n);

        % Plug into cell (in all fields)
        for k=1:length(g_file)
            if iscell(g_file{k})
                % insert empty string for cells
                if n == 1
                    % handling issue with the header
                    g_file{k} = {val_st, g_file{k}{:}};
                else
                    g_file{k} = {g_file{k}{1:idx},val_st,g_file{k}{idx+1:end}};
                end
            % for anything that isn't a cell (i.e. anything but the
            % dates/times)
            else
                % for the remaining columns (doubles):
                % fill with not a number (nan)
                if n == 1
                    % transpose so we avoid horzcat error
                    g_file{k} = [val_n,g_file{k}(:)'];
                else
                    g_file{k} = [g_file{k}(1:idx),val_n,g_file{k}(idx+1:end)];
                end
            end
        end
    end
    %----------------------------------------------------------------------
    % STORE DATES, TIMES, EVIDs, & MAGS FOR EACH EVENT IN CURRENT GRID
    %----------------------------------------------------------------------
    % variables to reset for each grid
    fam_number = 0; % family counter
   
    % create cells for each family in each grid
    fam_g = cell(length(g_empty),1); % empty cell to fill with values (time strings & evIDs)
    % n=1:length(g_file{1})
    for n=1:length(g_file{1})
        date_ev = g_file{1}{n}; % date for d'th cell
        time_ev = g_file{2}{n}; % time "
        id_ev = g_file{7}(n); % cluster id
        depth_ev = g_file{5}(n); % depth
        mag_ev = g_file{6}(n); % event magnitude
        % convert from coda duration magnitude to seismic moment (M0)
        m0_ev = 10^((1.2*mag_ev)+17); % (Bakun 1984) dyne-cm

        % check if date cell is empty, if yes, then assign number to family
             % will later have to remove families with recurrence intervals
             % <0.1 yrs
        if isempty(date_ev)
            % returns 1 if empty, 0 if not
            fam_number = fam_number+1; % assign family number
        else
            % if not empty, i.e. we have dates and times
            % append date and time for each event, d
            date_str = append(date_ev,' ',time_ev);
            date_time = datetime(date_str,'InputFormat','yyyy-MM-dd HH:mm:ss.SSSS');% stores fractional seconds but will not display them
            date_time.Format = 'yyyy-MM-dd HH:mm:ss.SSSS'; % to display fractional seconds
            
            % add to fam_g cell array
            fam_g{fam_number,1}{n,1} = date_time; % date time value
            % store magnitude
            fam_g{fam_number,1}{n,2} = m0_ev; % dyne-cm
            % store depth
            fam_g{fam_number,1}{n,3} = depth_ev; % km
            % store clusid
            fam_g{fam_number,1}{n,4} = id_ev;
            % end empty cell test; i.e. family counting    
        end
    end
    % get rid of empty cells
    for f = 1:length(fam_g)
        % remove empty times evIDs
        no_empties = fam_g{f}(~cellfun('isempty',fam_g{f})); % returns column
        
        % reshape no_empties to og size
        no_empties = reshape(no_empties,[(length(no_empties)/4),4]);
        % replace fam_g{f} with no_empties
        fam_g{f} = no_empties;
    end

    %***********************
    % Spring cleaning
    %***********************
    clear f n

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CALCULATE RECURRENCE INTERVALS (Tr)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculating recurrence intervals between consecutive events in a family
    %----------------------------------------------------------------------
    % Minimum recurrence interval
    % Chen (2008) excludes recurrence intervals <0.1 yrs from creep rate
    % estimates. See pg. 268; discusses possibility of these not being actual
    % repeaters, but part of triggering processes.
    % This might be an issue; after discussing with GF: RE Trs might speed
    % up following a bigger earthquake. Current code doesn't check for
    % that -- something to follow up on.
    %----------------------------------------------------------------------
    % CALCULATE RECURRENCE INTERVALS BETWEEN EACH EVENT IN A FAMILY
    % From Chen et al (2007):
    % If we have N events in a sequence ([i1 i2 i3 i4 ... iN]) then we will
    % have N-1 recurrence intervals, j ([(i2-i1) (i3-i2) (i4-i3)... (iN-i(N-1))])
    %----------------------------------------------------------------------
    % FILTER OUT CONSECUTIVE EVENTS WITH RECCURRENCE INTEVRALS <0.1 YRS
    fam_Tr = cell(length(fam_g),1); % new cell array to store updated list of families + their Tr
                                    % to be used in slip calculations (i.e. no Tr<Tr_min)
    fam_Tr_all = fam_Tr; % same as fam_Tr but it keeps ALL recurrence intervals
                         % to be used for alpha and beta calibrations or other analyses--
                         % TENTATIVE
    for f = 1:(length(fam_g))
        % sort events by date time
        % (some fams have events out of order -- issue from RE catalog)
        [~, idx] = sort([fam_g{f}{:,1}]); % sort by datetime 
        fam_g{f} = fam_g{f}(idx,:); % replace cell with sorted values
        
        % CALCULATE RECCURENCE INTERVAL (TR)
        t_fam = fam_g{f}(:,1); % times in family
        % Following eqn. 4 in Chen 2008:
        % Sequence has N events; i goes from 2 to N; j is recurrence
        % interval which goes from 1 to N-1
        for t=2:length(t_fam)
            t2 = t_fam{t}; % event i
            t1 = t_fam{t-1}; % event i-1
    
            t_diff = t2-t1; % recurrence interval; duration HH:MM:SS
            t_diff = seconds(duration(t_diff,'Format','hh:mm:ss.SSSS')); %convert to seconds
    
            % If Tr is >= 0.1 yrs, then store this info, along with M0,
            % depth, & clusId in the cell array to be used for slip calculations (fam_Tr)
            if t_diff >= Tr_min
                % slip cell array
                fam_Tr{f,1}{t,1} = t_diff; % t2-t1
                fam_Tr{f,1}{t,2} = [fam_g{f}{t,2},fam_g{f}{t-1,2}]; % M0 for both events (M0_2, M0_1)
                fam_Tr{f,1}{t,3} = [fam_g{f}{t,3},fam_g{f}{t-1,3}]; % depth for both events (D2, D2)
                fam_Tr{f,1}{t,4} = [fam_g{f}{t,4},fam_g{f}{t-1,4}]; % clusID for both events (clusID2, clusID1)
            end
    
            % store all recurrence intervals for calibrations in respective
            % cell array:
            % calibrations cell array
            fam_Tr_all{f,1}{t,1} = t_diff; % t2-t1; [sec]
            fam_Tr_all{f,1}{t,2} = [fam_g{f}{t,2},fam_g{f}{t-1,2}]; % M0 for both events (M0_2, M0_1) [dyne-cm] 
            fam_Tr_all{f,1}{t,3} = [fam_g{f}{t,3},fam_g{f}{t-1,3}]; % depth for both events (D2, D2) [km]
            fam_Tr_all{f,1}{t,4} = [fam_g{f}{t,4},fam_g{f}{t-1,4}]; % clusID for both events (clusID2, clusID1)
        end
    end
    %------------------------------------------------------
    % Getting rid of empty cells in fam_Tr and fam_Tr_all
    % this is VILE but it works.
    % TR > 0.1 YRS
    % Get rid of empty families first
    blah = {}; % temporary variable so I don't get confused
    b = 1;
    for f1=1:length(fam_Tr)
        if ~isempty(fam_Tr{f1})
            blah{b,1} = fam_Tr{f1};
            b=b+1;
        end
    end
    
    % get rid of empty Tr's within families
    for f1=1:length(blah)
        % leave non-empty cells
        no_empties = blah{f1}(~cellfun('isempty',blah{f1}));
        % reshape no_empties to og size
        no_empties = reshape(no_empties,[(length(no_empties)/4),4]);
        % replace fam{f} with no_empties
        blah{f1} = no_empties;
    end
    fam_Tr = blah;
    % ALL TR
    % for Tr all (i.e. Time difference between event i=1 and itself isn't
    % taken.)
    blah2 = fam_Tr_all; % for tr all; temporary variable so I don't get confused
    for f1=1:length(blah2)
        % leave non-empty cells
        no_empties = blah2{f1}(~cellfun('isempty',blah2{f1}));
        % reshape no_empties to og size
        no_empties = reshape(no_empties,[(length(no_empties)/4),4]);
        % replace fam{f} with no_empties
        blah2{f1} = no_empties;
    end
    fam_Tr_all=blah2;
    %------------------------------------------------------
    %***********************
    % spring cleaning
    %***********************
    clear b blah blah2 date_ev date_str date_time depth_ev f f1 id_ev...
        m0_ev mag_ev no_empties t t1 t2 t_diff T_sort time_ev

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% AVERAGE RECURRENCE INTERVAL & SEISMIC MOMENT FOR EACH SEQUENCE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fam_avg = cell(length(fam_Tr),1); % fam avg Tr and avg M0 (for slip calculations)
    fam_avg_all = cell(length(fam_Tr_all),1); % fam avg Tr and avg M0 (for calibrations)

    % for filtered Tr cell array:
    for f = 1:length(fam_Tr)
        % average Tr for family
        fam_avg{f,1} = mean([fam_Tr{f}{:,1}]);

        % average M0 for family
        M0s_fam = [fam_Tr{f}{:,2}];
        M0s_fam = [M0s_fam(2:2:end),M0s_fam(end-1)];
        fam_avg{f,2} = mean(M0s_fam);

        % average depth for family
        deps_fam = [fam_Tr{f}{:,3}];
        deps_fam = [deps_fam(2:2:end+1),deps_fam(end-1)];
        fam_avg{f,3} = mean(deps_fam); 
    end

    % for UNfiltered Tr cell array:
    for f = 1:length(fam_Tr_all)
        % average Tr for family
        fam_avg_all{f,1} = mean([fam_Tr_all{f}{:,1}]);

        % average M0 for family
        M0s_fam = [fam_Tr_all{f}{:,2}];
        M0s_fam = [M0s_fam(2:2:end),M0s_fam(end-1)];
        fam_avg_all{f,2} = mean(M0s_fam);

        % average depth for family
        deps_fam = [fam_Tr_all{f}{:,3}];
        deps_fam = [deps_fam(2:2:end+1),deps_fam(end-1)];
        fam_avg_all{f,3} = mean(deps_fam); 
    end

    % convert to mat
    fam_avg = cell2mat(fam_avg);
    fam_avg_all = cell2mat(fam_avg_all);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SAVE VARIABLES TO FILE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g_workspace = ['g_',num2str(g),'vars.mat'];
    save(g_workspace)

    %***********************
    % spring cleaning
    %***********************
    clear M0s_fam deps_fam empty_str f fam_avg fam_avg_all fam_g fam_number...
        fam_Tr fam_Tr_all fileID g g_empty g_file g_workspace idx k ...
        t_fam
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD TR, M0 FOR ALL GRIDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beta_calib takes in Tr and M0
% Limiting to Tr >= 0.1 yrs
% (Chen 2008) : AVERAGE Tr value and M0 value for each sequence

% ARRAYS TO APPEND TR AND M0 VALUES TO
Tr = []; % avg. reccurrence interval (sec) for ea. family on all grids (Tr >=0.1 yrs)
M0 = []; % avg. seismic moment (dyne-cm) for ea. family on all grids (Tr >=0.1 yrs)

Tr_all = []; % reccurrence interval (sec) for families on all grids (Tr all)
M0_all = []; % seismic moment (dyne-cm) for families on all grids (Tr all)

% loop through grids
for g=1:num_gs
    % load workspace
    g_workspace = ['g_',num2str(g),'vars.mat'];
    load(g_workspace)
    
    % loop through families in grid
    for f=1:length(fam_Tr)
        Tr = [Tr, fam_avg(f,1)]; % sec
        M0 = [M0, fam_avg(f,2)]; % dyne-cm
    end

    for f=1:length(fam_Tr_all)
        Tr_all = [Tr_all, fam_avg_all(f,1)]; % sec 
        M0_all = [M0_all, fam_avg_all(f,2)]; % dyne-cm
    end
end

%***********************
% spring cleaning
%***********************
clear M0s_fam deps_fam empty_str f fam_avg fam_avg_all fam_g fam_number...
    fam_Tr fam_Tr_all fileID g g_empty g_file g_workspace idx k ...
    t_fam
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BETA CALIBRATION + JACKKNIFE TEST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nadeau (1998): log(Tr) = alpha_T + beta_T*log(M0) (eqn. 15)
% Tr: avg. recurrence interval of a sequence
% If the average rate of seismic slip is indpendent of M0:
%       beta_T = beta_d  (Nadeau 1998, eqn. 19)
% If this is not the case, then eqn. 15 is an approximation.
%--------------------------------------------------------------------------
% BETA SLIP PARAMETER CALIBRATION
%--------------------------------------------------------------------------
% Plug Tr and M0 into alpha_calib(Tr_in,M0_in) for Tr in seconds:
[x] = beta_calib([Tr',M0']); % one alpha value for the entire fault
beta_T = x(1); alpha_T = x(2); % beta_T = beta_d
%--------------------------------------------------------------------------
% JACKKNIFE TEST (BETA_T)
% Removes 1 sequence (family) average value from list of all families (Tr>=0.1 yrs).
%--------------------------------------------------------------------------
%jackstat = jackknife(jackfun,x,y...)
% jackfun = function to perform on X
% From 'help':
% "JACKSTAT = jackknife(JACKFUN,X) draws jackknife data samples from the
%     N-by-P data array X, computes statistics on each sample using the
%     function JACKFUN, and returns the results in the matrix JACKSTAT.
%     JACKFUN is a function handle specified with @.  Each of the N rows of
%     JACKSTAT contains the results of applying JACKFUN to one jackknife
%     sample.  Row I of JACKSTAT contains the results for the sample
%     consisting of X with the Ith row omitted:"
%--------------------------------------------------------------------------
x_jack = jackknife(@beta_calib,[Tr',M0']); % array of beta and b_const calibrations for each test
                                            % col1 = beta values
                                            % col2 = b_const. values
beta_jack = x_jack(:,1); alpha_T_jack = x_jack(:,2); % beta_T = beta_d

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT JACKKNIFE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig_jack = figure(1);
scatter([1:length(Tr)],beta_jack,'b.','DisplayName','jackknife')
hold on
yline(beta_T,'LineWidth',2.5,'DisplayName',['\beta_{LSQ} = ',num2str(beta_T)]) % beta value from least-squares
xlabel('family removed (index)')
ylabel('\beta_T value')
title('\beta jackknife sampling results')
subtitle(['mean \beta=',num2str(mean(beta_jack)),'; std = ',num2str(std(beta_jack))])
legend()
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT LOG10(TR) VS LOG10(M0) AND BEST FIT LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try to match Fig. 11 in Nadeau (1998)
% avg Tr. and M0 for ea. family on the entire fault polygon
fig_logTRM0 = figure(2);
% plot all events
p1 = scatter(log10(M0_all),log10(Tr_all),5,'MarkerEdgecolor',[0.3 0.3 0.3],...
    'MarkerFaceColor','none');
hold on
% plot events with Tr>=0.1yrs
p2 = scatter(log10(M0),log10(Tr),5,'MarkerEdgecolor','b',...
    'MarkerFaceColor','b');

% beta calibration
x_vals = linspace(min(M0),max(M0)*10,500);
y_vals = (10^alpha_T)*(x_vals.^beta_T); % Nadeau (1998) eqn. 15: log(Tr)=alpha_T+beta_T*log(M0)
                                         % --> Tr = (10^alpha_T)*(M0^beta_T)
p3 = plot(log10(x_vals),log10(y_vals),'k','LineWidth',1.5);

% Parkfield beta (Nadeau 1998, eqn. 16)
y_park = (10^4.85)*(x_vals.^betaP);
p4 = plot(log10(x_vals),log10(y_park),'k','LineStyle','--');

% plot accoutrements
xlim([17 21])
ylim([-1 10])
title('\beta calibration  for Tr_{avg} \geq 0.1 yrs (Hayward Fault)')
xlabel('log M_0 (dyne-cm)')
ylabel('log Tr (sec)')
legend({['Tr_{avg} (all)'],['Tr_{avg} \geq 0.1yrs'] ...
    ['\alpha_{LSQ} = ' num2str(alpha_T) newline '\beta_{LSQ} = ' num2str(num2str(beta_T))],...
    ['\alpha_{Nadeau98} = ' num2str(4.85) newline '\beta_{Nadeau98} = ' num2str(betaP)]})
legend('Location','southeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USGS CREEP RATE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From: https://www.usgs.gov/data/creep-rate-models-california-faults-2023-us-national-seismic-hazard-model
%
% Johnson, K., Murray, J.R., and Wespestad, C.E., 2022, 
% Creep rate models for California faults in the 2023 US National Seismic 
% Hazard Model: U.S. Geological Survey data release, https://doi.org/10.5066/P94YGVWQ.
%--------------------------------------------------------------------------
% "The updated surface creep rate compilation consists of variety of data 
% types including alignment arrays, offset cultural markers, creepmeters, 
% InSAR, and GPS data. We compile a total of 497 surface creep rate 
% measurements, 400 of which are new and 97 of which appear in the UCERF3 
% compilation."
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD CREEP RATE MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
creep_models = readtable('creep_models.csv');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VALUES FOR HAYWARD FAULT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx = find(strcmp('Hayward',creep_models.Fault)==1); % find row indices of Hayward Fault values

lon = creep_models.Longitude(idx); % longitude
lat = creep_models.Latitude(idx); % latitude
yrs = creep_models.Years(idx); % total years; some nan values in here!

creep_rate = creep_models.CreepRate_mm_yr_(idx); % creep rate (mm/yr)
creep_rate = creep_rate/10; % mm to cm
creep_rate = creep_rate/(60*60*24*365); % cm/yr --> cm/sec
avg_creep = mean(creep_rate); % cm/sec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE COORDS. TO PLOT ON GOOGLE EARTH ON TOP OF HAYWARD FAULT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
coords = [lon,lat];
writematrix(coords,'HF_creep_model.txt')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STATIONS INSIDE HAYWARD POLYGONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find creep measurents located inside Hayward polygon
grids = load("grids.txt");
g_lons = grids(1:2:end,:);
g_lats = grids(2:2:end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GRID 13 CREEP RATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% from gmt plot (produced by plot_creep_model.gmt6) we can see the creep
% rate measurement locations from the USGS file and the RE failies are very
% well aligned (for the most part).
%--------------------------------------------------------------------------
% [in,on] = inpolygon(xq,yq,xv,yv) also returns on indicating if the query 
%                                   points are on the edge of the polygon area.
% in = inpolygon(xq,yq,xv,yv) returns in indicating if the query points 
%                             specified by xq and yq are inside or on the 
%                             edge of the polygon area defined by xv and yv.
%--------------------------------------------------------------------------
[in, on] = inpolygon(lon,lat,g_lons(13,:),g_lats(13,:)); % creep model rates inside grid 13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ALPHA_D ESTIMATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do I feel good about this part? Not in the least. But I'm running very
% short on time. SO:
% Nadeau 1998:
% eqn. 17: d_avg_slip_rate = d/Tr
% eqn. 20: alpha_T = alpha_d - log(d_avg_slip_rate)
% rearrange: alpha_d = alpha_T + log(d_avg_slip_rate)
%--------------------------------------------------------------------------
alpha_d = alpha_T + log10(mean(creep_rate(in==1)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SLIP AND CREEP CALCULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g_d_avgs = cell(num_gs,1); % cell array containing avg slips for each family in a grid
g_M0_avgs = cell(num_gs,1); % cell array containing avg M0 for each family in a grid
g_Vd_short = cell(num_gs,1); % array containing avg Vd_short for ea. family in a grid
g_Vd_long = cell(num_gs,1); % array containing avg Vd_long for ea. family in a grid

for g = 1:num_gs
    % ALPHA & BETA FROM THIS PROJECT
    % CALCULATE SLIP & CREEP RATES
    [di_mean, M0_mean, Vd_short, Vd_long] = slip_creep(g,alpha_d,beta_T);

    % SAVE VALUES - FROM THIS PROJECT
    g_d_avgs{g,1} = di_mean;
    g_M0_avgs{g,1} = M0_mean;
    g_Vd_short{g,1} = Vd_short; 
    g_Vd_long{g,1} = Vd_long;

    % ALPHA & BETA FROM PARKFIELD (NADEAU 1998)
    % CALCULATE SLIP & CREEP RATES
    [di_mean, M0_mean, Vd_short, Vd_long] = slip_creep(g,alphaP,betaP);

    % SAVE VALUES - FROM NADEAU98
    g_d_avgs{g,2} = di_mean;
    g_M0_avgs{g,2} = M0_mean;
    g_Vd_short{g,2} = Vd_short; 
    g_Vd_long{g,2} = Vd_long;

    % ALPHA & BETA FROM KHOSHMANESH (2015)
    % CALCULATE SLIP & CREEP RATES
    [di_mean, M0_mean, Vd_short, Vd_long] = slip_creep(g,alphaK,betaK);

    % SAVE VALUES - FROM KHOSHMANESH (2015)
    g_d_avgs{g,3} = di_mean;
    g_M0_avgs{g,3} = M0_mean;
    g_Vd_short{g,3} = Vd_short; 
    g_Vd_long{g,3} = Vd_long;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT SLIP VS SEISMIC MOMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting slip vs seismic moment for all the sequences in the fault
% Comparing these results to Nadeau 1998's alpha and beta slip.
%--------------------------------------------------------------------------
fig_D_M0 = figure(3);
sgtitle(['d_i vs. M_0 comparison'])
for g=1:num_gs
    subplot(3,5,g); % create subplot
    % this project's results
    p1 = scatter(log10(g_M0_avgs{g,1}),log10(g_d_avgs{g,1}),'b','.');
    hold on

    % Parkfield (Nadeau 1998)
    p2 = scatter(log10(g_M0_avgs{g,2}),log10(g_d_avgs{g,2}),'r','.');

    % Khoshanesh 2015
    p3 = scatter(log10(g_M0_avgs{g,3}),log10(g_d_avgs{g,3}),'m','.');

    % plot stuff
    xlabel('log M_{0} (dyne-cm)')
    ylabel('log d (cm)')
    title(['Grid ',num2str(g)])
end
% add overall legend
lh = subplot(3,5,15);
set(get(lh,'Children'),'Visible','off');
% fake points just to pot soemthing for the legend
plot(lh,1,1,'b.');
hold on
plot(lh,1,1,'r.');
plot(lh,1,1,'m.');
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
legend({['\alpha_d = ' num2str(alpha_d) newline '\beta_d = ' num2str(num2str(beta_T))],...
    ['\alpha_{Nadeau98} = ' num2str(alphaP) newline '\beta_{Nadeau98} = ' num2str(betaP)]...
    ['\alpha_{Khosh15} = ' num2str(alphaK) newline '\beta_{Khosh15} = ' num2str(betaK)]})
legend('Location','north')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHORT-TERM CREEP RATE COMPARISONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USGS AVG. SHORT-TERM CREEP RATES NOT AVAILABLE
g_Vd_short_avg = zeros(num_gs,3); % avg. Vd long for each grid (all 3 calibrations)

for g=1:num_gs
    % find avg. Vd_short for each grid
    g_Vd_short_avg(g,1) = mean(g_Vd_short{g,1}); % mine
    g_Vd_short_avg(g,2) = mean(g_Vd_short{g,2}); % Nadeau
    g_Vd_short_avg(g,3) = mean(g_Vd_short{g,3}); % Khosh
end
% convert from cm/sec to cm/yr
% mine & nadeau
g_Vd_short_avg = g_Vd_short_avg*(60*60*24*365);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT VD_SHORT COMPARISONS (CM/YR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = [1:num_gs];
fig_Vdshort = figure(4);
% mine
scatter(g,g_Vd_short_avg(:,1),'k','filled','DisplayName',['\alpha_d=',num2str(alpha_d),'; \beta_d=',num2str(beta_T)])
hold on
% Nadeau
scatter(g,g_Vd_short_avg(:,2),'k','LineWidth',1.2,...
    'DisplayName',['\alpha _{Nad98} =',num2str(alphaP),'; \beta _{Nad98} =',num2str(betaP)])
% Khoshmanesh
scatter(g,g_Vd_short_avg(:,3),'k','s','LineWidth',1.2,...
    'DisplayName',['\alpha _{Khosh15} =',num2str(alphaK),'; \beta _{Khosh15} =',num2str(betaK)])

xlabel('Grid number')
ylabel("V_d (cm/yr)")
title('Short-term creep rate comparison')
legend()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LONG-TERM CREEP RATE COMPARISONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USGS AVG. LONG-TERM CREEP RATE
in_poly = zeros(num_gs,1);
on_poly = zeros(num_gs,1);
g_Vd_long_avg = zeros(num_gs,3); % avg. Vd long for each grid (all calibrations)

for g=1:num_gs
    %----------------------------------------------------------------------
    % find creep rates inside grids
    [in on] = inpolygon(lon,lat,g_lons(g,:),g_lats(g,:));

    if length(find(in==1))>0
        % not all grids have rates available, so check for that.
        in_poly(g) = mean(creep_rate(find(in==1))); % cm/yr
    else length(find(on==1))>0;
        on_poly(g) = mean(creep_rate(find(on==1))); % cm/yr
    end
    %----------------------------------------------------------------------
    % find avg. Vd_long for each grid while we're looping through grids
    g_Vd_long_avg(g,1) = mean(g_Vd_long{g,1}); % mine
    g_Vd_long_avg(g,2) = mean(g_Vd_long{g,2}); % Nadeau
    g_Vd_long_avg(g,3) = mean(g_Vd_long{g,3}); % Khoshmanesh
end
% convert from cm/sec to cm/yr
% usgs
in_poly = in_poly*(60*60*24*365);
% mine & nadeau
g_Vd_long_avg = g_Vd_long_avg*(60*60*24*365);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT VD_LONG COMPARISONS (CM/YR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = [1:num_gs];
fig_VdLong = figure(5);
% mine
scatter(g,g_Vd_long_avg(:,1),'k','filled','DisplayName',['\alpha_d=',num2str(alpha_d),'; \beta_d=',num2str(beta_T)])
hold on
% Nadeau
scatter(g,g_Vd_long_avg(:,2),'k','LineWidth',1.2,...
    'DisplayName',['\alpha _{Nad98} =',num2str(alphaP),'; \beta _{Nad98} =',num2str(betaP)])
% Khoshmanesh
scatter(g,g_Vd_long_avg(:,3),'k','s','LineWidth',1.2,...
    'DisplayName',['\alpha _{Khosh15} =',num2str(alphaK),'; \beta _{Khosh15} =',num2str(betaK)])
% USGS creep model - only found in grids 10-14
scatter(g(10:14),in_poly(10:14),40,'k','^','filled','DisplayName','USGS creep model')

xlabel('Grid number')
ylabel("V_d (cm/yr)")
title('Long-term creep rate comparison')
legend()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPORT SHORT-TERM CREEP RATES FOR MY VALUES
writematrix(g_Vd_short_avg(:,1),'short_term_creep.txt')
% EXPORT LONG-TERM CREEP RATES FOR MY VALUES
writematrix(g_Vd_long_avg(:,1),'long_term_creep.txt')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE FIGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% but first: set figure sizes
% https://www.mathworks.com/matlabcentral/answers/65402-how-to-set-graph-size#answer_76937
llc_x=20; % lower left corner (llc) of screen to llc of figure
llc_y=20; % lower left corner (llc) of screen to llc of figure
w=800; % figure width
h=700; % figure height

% fig1
set(fig_jack,'position',[llc_x,llc_y,w,h])
saveas(fig_jack,'jackknife_results.png')
% fig2
set(fig_logTRM0,'position',[llc_x,llc_y,w,h])
saveas(fig_logTRM0,'Tr_M0_comparison.png')
% fig3
set(fig_D_M0,'position',[llc_x,llc_y,1200,h])
saveas(fig_D_M0,'slip_M0_comparison.png')
%fig4
set(fig_Vdshort,'position',[llc_x,llc_y,w,400])
saveas(fig_Vdshort,'short_term_creep.png')
%fig5
set(fig_VdLong,'position',[llc_x,llc_y,w,400])
saveas(fig_VdLong,'long_term_creep.png')