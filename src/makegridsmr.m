
function [pdb,smras,smrvs,gsta,stinfo,VphaseMod]=makegridsmr(par,isfigure,ddst,remakegrd,prdb,smrv,smra)
%% ------- %% make smra smrv--------
pbeg=min(par.periods)-1;
pend=max(par.periods)+1;
pdb=[pbeg prdb pend]; %period boundary for smoothing radius
smrvs(length(pdb)-1)=0;
smras(length(pdb)-1)=0; 
for i=1:length(pdb)-1
    id1=par.periods>pdb(i)&par.periods<=pdb(i+1)+1;
    smrvs(id1)=smrv(i);
    smras(id1)=smra(i);
end
% plot(par.periods,smrvs)
% hold on
% plot(par.periods,smras)
%% ------read station information
[~,stnm,stloc(:,1),stloc(:,2),stloc(:,3)]=textread(par.stinfof,'%s%s%f%f%f'); %#ok<DTXTRD>
stinfo.stnm = stnm;
stinfo.stloc = stloc;
%% ----- make or read vphase mod
if isempty(par.vmodtype)
VphaseMod=load(par.v_modf);
VphaseMod=VphaseMod.VphaseMod;
elseif strcmp(par.vmodtype,'ak135')
T=par.periods;    
freq=1./T;
freq= sort(freq,'ascend');
vphase=ak135_dispersion(freq);
vphase=sort(vphase,'ascend');
VphaseMod=[vphase',T'];
end

%% ------- make grid file ----------
if ~exist('gstall.mat','file')||remakegrd==1
    k=0;    
    for i=par.lalim(1):par.gridsize:par.lalim(2)
        for j=par.lolim(1):par.gridsize:par.lolim(2)
            lats=stloc(stloc(:,1)<i+ddst & stloc(:,1)>i-ddst &stloc(:,2)<j+ddst & stloc(:,2)>j-ddst,1);
            lons=stloc(stloc(:,1)<i+ddst & stloc(:,1)>i-ddst &stloc(:,2)<j+ddst & stloc(:,2)>j-ddst,2);
            evls=stloc(stloc(:,1)<i+ddst & stloc(:,1)>i-ddst &stloc(:,2)<j+ddst & stloc(:,2)>j-ddst,3);
            if length(lats)<par.minNST
                continue
            end
            dist=distance(lats,lons,i,j);
            minID = dist==min(dist);            
            gevl=mean(evls(minID));            
            k=k+1;
            st(k,:)=[round(i*100)/100,round(j*100)/100,round(gevl*100)/100]; %#ok<*AGROW,*SAGROW>
            gsta(k).st=st(k,:);
            gsta(k).stn=['bk' num2str2(k,4,0)];
        end
    end
    save('gstall.mat','gsta')
else
    load('gstall.mat')
    st=cell2mat({gsta(:).st}'); %#ok<*NODEF>
end

if isfigure==1
    figure(595)
    plot(stloc(:,2),stloc(:,1),'^k')
    hold on
    plot(st(:,2),st(:,1),'.r')
    axis equal
    hold off
    drawnow
end