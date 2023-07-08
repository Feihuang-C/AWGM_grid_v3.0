%--------------------------------------------
%---irregular array----
%--------------------------------------------
function [dudx,dudy,stc]=fitVg(st,ev,fband,par,sc1)
vavg=par.vavg;
damp=0.01;
fc=mean(fband); 
nst=length(st);
    %---filter the waves---
np=length(st(sc1).dat);
    %Form the equation system to solve for spatial gradient
    %   d=Gx
stc=st(1);
dt=st(2).dt;

for is=2:nst
    [envu,~]=envelope(st(is).dat);
    kmax=findmaxima(envu,1);
    ktime(is-1)=kmax*dt+st(is).tbeg;
end

k=0;
for is=2:nst
        ds1=distance(st(is).st(1),st(is).st(2),st(is).ev(1),st(is).ev(2));
        azm1=azimuth(st(is).st(1),st(is).st(2),st(is).ev(1),st(is).ev(2));
    for js=2:nst
        k=k+1;
        azm2=azimuth(st(is).st(1),st(is).st(2),st(js).st(1),st(js).st(2));
        ds2=distance(st(js).st(1),st(js).st(2),st(js).ev(1),st(js).ev(2));
        dds=ds1-ds2;
        ds=distance(st(is).st(1),st(is).st(2),st(js).st(1),st(js).st(2))*111.1949*cos((azm2-azm1)/180*pi)*dds/abs(dds);
        wi(k)=compweight(azm1,st(1),st(is),damp,par.kwi);
        dt=ktime(is-1)-ktime(js-1); 
        vgsum(k)=ds/dt;
    end
end

vg=sum(vgsum(~isnan(vgsum)))/sum(wi(~isnan(vgsum)));

nr=length(dx);
if par.kwi>0
    for is=1:nr
        wii=wi(is);
        if ~isempty(wii)
            ui(is,1:np)=ui(is,1:np)*wii;
            dx(is)=dx(is)*wii;
            dy(is)=dy(is)*wii;
            
        end
    end
end
G=[wi' dx' dy'];
u0=(pinv(G'*G))*G'*ui;
stc(1).dat=u0(1,:);
dudx=u0(2,:);
dudy=u0(3,:);


function [kmax]=findmaxima(dat,mmax)
diff0=diff(dat);
ndf=length(diff0);
nmax=0;
for id=2:ndf-1
    if diff0(id)==0 & diff0(id-1)>0 & diff0(id+1)<0
        nmax=nmax+1;
        pmax(nmax)=id;
        %[blw(nmax),bup(nmax)]=errbars(dat,id,0.95);
    elseif diff0(id)>0 & diff0(id+1)<0
        nmax=nmax+1;
        pmax(nmax)=id+1;
        %[blw(nmax),bup(nmax)]=errbars(dat,id+1,0.95);
    end
end

if nmax==0
    kmax=[];
    return
end

mmax=min(length(pmax),mmax);
for ii=1:mmax
    [fmax,kk]=max(dat(pmax));
    kmax(ii)=pmax(kk);
    pmax(kk)=[];
end

%%
%% Calculate the weights for station pairs.
%%
function wi=compweight(azpth,st1,st2,damp,kwi)
if ~isempty(azpth)
    azst=azimuth0(st1.st(1),st1.st(2),st2.st(1),st2.st(2));
    dst=distance0(st1.st(1),st1.st(2),st2.st(1),st2.st(2))*111.149;
    daz=(azpth-azst)/180*pi;
    wi=dst*cos(daz);
    wi=0.1/(abs(wi)+damp);
    wi=wi^kwi;
    %wi=wi*abs(sin(daz)^(kwi));
else
    wi=[];
end

%%
%% Calculate the errors for truncation of first order taylor series.
%%
function delta=CompErrors(azpth,st1,st2,vel,fc)
if ~isempty(azpth)
    azst=azimuth0(st1.st(1),st1.st(2),st2.st(1),st2.st(2));
    dst=distance0(st1.st(1),st1.st(2),st2.st(1),st2.st(2)); %in degree
    daz=azpth-azst;
    wi=dst*cos(daz/180.0*pi);
    lamda=vel/fc/111.149;
    delta=pi/lamda*abs(wi);
else
    delta=[];
end

function [CS,isbad] = CSestimate(u0,st,fband,TimeMaxVg,par)
    st(1).dat = u0;
    periods = sort(1./fband);
    period = mean(periods);   
    % Calculate Auto-correlation
    ijk = 0;
for ista = 1:length(st)
    stlas(ista)=st(ista).st(1);
    stlos(ista)=st(ista).st(2);
    autocor = CSestimate0(st,ista,ista,periods,TimeMaxVg,period,par.dt);
    if sum(autocor.amp) == 0
        ijk = ijk+1;
        badID(ijk) = ista;
    else
        CSx.autocor(ista) = autocor;
        badID = [];
    end
end
if ~isempty(badID)
st(badID) = [];

if length(st)<4 ||find(ismember(badID,1))>0
    CS = [];
    return
end
end
% start to find nearby stations and apply cross-correlation measurement
for ista = 2:length(st)
     CS(ista-1) = CSestimate0(st,1,ista,periods,TimeMaxVg,period,par.dt);
end % end of station loop


if ~exist('CSx','var')
    CS = [];
end

% Calculate the coherency between each station pairs
for ics = 1:length(CS)
    sta1 = CS(ics).sta1;
    sta2 = CS(ics).sta2;
    CS(ics).cohere = CS(ics).amp^2/CSx.autocor(sta1).amp/CSx.autocor(sta2).amp;
end

for ics = 1:length(CS)
    sta1 = CS(ics).sta1;
    sta2 = CS(ics).sta2;
    CS(ics).cohere = CS(ics).amp^2/CSx.autocor(sta1).amp/CSx.autocor(sta2).amp;
end

for ics = 1:length(CS)
    %remove the station pairs with low coherency
    if CS(ics).cohere>=0.5
        CS(ics).isgood = 1;
    else
        CS(ics).isgood = -1;
    end
    %remove the station pairs with great error in dtp
%     if abs(CS(ics).dtp-median(detrend([CS.dtp])))<=10
%         CS(ics).isgood = 1;
%     else
%         CS(ics).isgood = -1;
%     end
end
isbad = find([CS.isgood]<0)+1;
CS([CS.isgood]<0)=[];

if length(CS)<4 
    CS = [];
    return
end


%find the subarray for each station

    
function CS = CSestimate0(sta,sta1,sta2,periods,TimeMaxVg,cntprd,dt)
    CS.sta1 = sta1;
    CS.sta2 = sta2;
    CS.stan1 = sta(sta1).stn;
    CS.stan2 = sta(sta2).stn;


    data1 = sta(sta1).dat;
    bgtime = sta(sta1).tbeg;
    dt1 = sta(sta1).dt;
    Nt = length(data1);
    fN = 1/2/dt1;
%     [b,a] = butter(2,[1/(periods(end))/fN, 1/(periods(1))/fN]);
%     data1 = filtfilt(b,a,data1);
    taxis1 = bgtime + [0:Nt-1]'*dt1;
    dist1 = sta(sta1).dis;


    data2 = sta(sta2).dat;
    bgtime = sta(sta2).tbeg;
    dt2 = sta(sta2).dt;
    Nt = length(data2);
    fN = 1/2/dt2;
%     [b,a] = butter(2,[1/(periods(end))/fN, 1/(periods(1))/fN]);
%     data2 = filtfilt(b,a,data2);
    taxis2 = bgtime + [0:Nt-1]'*dt2;
    dist2 = sta(sta2).dis;

    % resample the data if necessary
if dt1 > dt2
	new_taxis2 = taxis2(1):dt1:taxis2(end);
	data2 = interp1(taxis2,data2,new_taxis2);
	taxis2 = new_taxis2;
	dt2 = dt1;
elseif dt1 < dt2
	new_taxis1 = taxis1(1):dt2:taxis1(end);
	data1 = interp1(taxis1,data1,new_taxis1);
	taxis1 = new_taxis1;
	dt1 = dt2;
end

% apply cross-correlation
[xcor,lag] = xcorr(data1,data2,...
             floor(length(data1)*dt));
lag = lag.*dt1;
lag = lag + taxis1(1) - taxis2(1);
[para,resnorm,residual, exitflag] = gsdffit(xcor,lag,1./cntprd,2);
CS.dtp = para(4);
CS.dtg = para(5);
CS.amp = para(1);
CS.w = para(2);
CS.sigma = para(3);
CS.dis = dist1-dist2;
