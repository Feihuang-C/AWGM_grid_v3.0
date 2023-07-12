%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function par=aa_WGM_parameters
clc
clear;
isfigure=1; 
%% ----- 1. info of seismic array ------------
par.fs = '/';
par.stinfof = './st.info'; % station information: Net Station Name stla stlo stevl
par.ev_sac_rootpath = './sacdata/'; % root path of sac files 
par.lalim = [26 32.4];  % lat range
par.lolim = [99.6 105.4]; % lon range
par.cmp = 'BHZ';    % component
par.Net='T1';       % NET
par.gridsize = 0.2;   % deg, grid size in degrees
ddst=0.7;             % deg, isolate the reference locations (grid points) with less supporting stations  
remakegrd = 1;        % =1, remake grid file

%% ----- 2. Auto waveform pick --------
par.wvPick_output = './wvPick/'; %output path 
par.periods = [20 40]; % s, central periods s
par.ev_minDist = deg2km(2);     % km, minimum epicentral distance
par.ev_maxDist = deg2km(160);   % km, maximum epicentral distance
par.ev_maxDepth = 100;  % km, maximum depth of earthquakes
par.vglen = [2,5];      % km/s, minimum and maximum Vg (km/s) to isolate Surface wave
par.dvg = 0.2;          % km/s,outlier of Surface wave group arrival
par.SNR=5;              % minimum SNR of one waveform
par.mSNR=8;             % minimum mean(SNR) of one earthquake

%% ------- 3. parameters of WGM --------
%               subarray parameters  
par.WG_output = './WG/'; % output path 
par.minNST=6;         % minimum number of stations in a subarray (NSS)
par.maxstadist = 30; % km, subarray maxR
par.maxdsmax = 150; % km, maximum maxR when NSS<miNST
par.cutWin = 500; % s, time window around the surface wave
par.dvmax = 0.01; % km/s, min dv in the reducing velocity method
par.Rcof=0.7; % minimum correlation coefficient between the predicted wavefield and obseved wavefield

%                smoothing parameters 
par.vmodtype=[ ];% ='ak135',reference phase velocity from AK135; =[],user difine reference phase velocity, 
par.v_modf='./VphaseMod.mat'; % when par.vmodtype=[], user needs to give reference-velocity file here
par.vlim=0.8; % outlier<median(v)-vlim  or  outlier > median(v)+vlim
prdb=[18 50 60];      		%period boundary for smoothing radius
smrv=[0.4 0.3 0.4 0.5]; %deg,smoothing radius of velocity
smra=[0.6 0.5 0.6 0.7]; %deg,smoothing radius of anisotropy
%% ---------------------------------------------
[pdb,smras,smrvs,gsta,stinfo,vphasemod]=makegridsmr(par,isfigure,ddst,remakegrd,prdb,smrv,smra);
par.pdb=pdb;
par.smras=smras;
par.smrvs=smrvs;
par.gsta=gsta;
par.stinfo=stinfo;
par.vmod=vphasemod;







