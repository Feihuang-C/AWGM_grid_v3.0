function varargout = ae_step4_AnisoFitGUI(varargin)
% AE_STEP4_ANISOFITGUI MATLAB code for ae_step4_AnisoFitGUI.fig
%      AE_STEP4_ANISOFITGUI, by itself, creates a new AE_STEP4_ANISOFITGUI or raises the existing
%      singleton*.
%
%      H = AE_STEP4_ANISOFITGUI returns the handle to a new AE_STEP4_ANISOFITGUI or the handle to
%      the existing singleton*.
%
%      AE_STEP4_ANISOFITGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AE_STEP4_ANISOFITGUI.M with the given input arguments.
%
%      AE_STEP4_ANISOFITGUI('Property','Value',...) creates a new AE_STEP4_ANISOFITGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ae_step4_AnisoFitGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ae_step4_AnisoFitGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ae_step4_AnisoFitGUI

% Last Modified by GUIDE v2.5 30-Jun-2023 08:44:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ae_step4_AnisoFitGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ae_step4_AnisoFitGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ae_step4_AnisoFitGUI is made visible.
function ae_step4_AnisoFitGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ae_step4_AnisoFitGUI (see VARARGIN)

% Choose default command line output for ae_step4_AnisoFitGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ae_step4_AnisoFitGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

par=getparameters(handles);
par.fs='/';
% stas=dir([par.indir par.fs '*.mat']);
% par.stas=stas;
pstr(1:length(par.periods))={nan};
k=0;
for i=par.periods
    k=k+1;
    pstr(k)={i};
end
par.period=par.periods(1);
set(handles.PeriodList,'string',pstr);
set(handles.PeriodList,'value',1);
% set(handles.stationList,'string',{stas.name});
try
FF=get(handles.faultf,'TooltipString');
FF=strsplit(FF,':');
par.FF=cell2mat(FF(2));

BBF=get(handles.BBF,'TooltipString');
BBF=strsplit(BBF,':');
par.BBF=cell2mat(BBF(2));

faults=readfault(par.FF);
Blocks=readfault(par.BBF);
par.faults=faults;
par.Blocks=Blocks;
subplot(handles.Plot5);
hold off
for bi=1:length(Blocks)
    if ~isempty(Blocks(bi).a)
        plot(Blocks(bi).a(:,1), Blocks(bi).a(:,2),'color',[0.1660 0.540 0.1880],'Linewidth', 2)
        hold on
    end
end

for bi=1:length(faults)
    if ~isempty(faults(bi).a)
        plot(faults(bi).a(:,1), faults(bi).a(:,2),'k','Linewidth', 0.5)
    end
end
xlim(par.Lonlim)
ylim(par.Latlim)
% axis equal
catch
end
userdata.par=par;
set(gcf,'userdata',userdata);


% --- Outputs from this function are returned to the command line.
function varargout = ae_step4_AnisoFitGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in stationList.
function stationList_Callback(hObject, eventdata, handles)
% hObject    handle to stationList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns stationList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from stationList
userdata=get(gcf,'userdata');
par=userdata.par;
par.fs='/';
fout0= get(hObject,'String');
par.ksta=fout0(get(hObject,'Value'),:);
par.cdir=[par.indir  par.fs cell2mat(par.ksta)];
vlim=par.vlim;
isCI=get(handles.CIR,'value');
azmbin=get(handles.azmbin,'string');
azmbin=str2double(azmbin);
wgst=load(par.cdir);
wgst=wgst.wgst;
evla=wgst.evla;
evlo=wgst.evlo;
st=wgst.st;
pdx=wgst.periods==par.period;
dazm=wgst.dazm(:,pdx);
azmo=wgst.azm(:,pdx);
azmo(azmo<0)=azmo(azmo<0)+360;
azmo(azmo>360)=azmo(azmo>360)-360;
va=wgst.va(:,pdx);
v=wgst.v(:,pdx);
vf=va;
azmf=azmo;

if isempty(va)
    disp('no data')
    return
end

% remove bad velue
ind= isnan(v) | isnan(va) |abs(dazm)>30 | v==0 | va==0;

% ind=isnan(v);
azmo(ind)=[];
va(ind)=[];
dazm(ind)=[];
evla(ind)=[];
evlo(ind)=[];
viso=median(va);

if isempty(va)
    disp('no data')
    return
end

% remove outlier 1
ind=abs(va-viso)>vlim;
azmo(ind)=[];
va(ind)=[];
dazm(ind)=[];
evla(ind)=[];
evlo(ind)=[];
viso=median(va);

if length(va)<5
    disp('no enough data,skip')
end

% remove outlier 2
if isCI==1
    An1 = AnisoFit(azmbin,va,azmo,360,0);
    vfit=fun1(An1.pout,azmo/180*pi);
    dv=va-vfit;
    ind=abs(dv)>vlim*0.5;
    azmo(ind)=[];
    va(ind)=[];
    dazm(ind)=[];
    evla(ind)=[];
    evlo(ind)=[];
    viso=median(va);
end
%% for test
% azm0=0:360;
% v=viso+0.03*cos(2*(azm0-45)/180*pi)+0.1*cos(azm0/180*pi-33);
% vf=v;
% azmf=azm0;
%%
azm2=azmo;
while ~isempty(find(azm2>180|azm2<0,1))
    azm2(azm2>180) = azm2(azm2>180)-180;
    azm2(azm2<0) = azm2(azm2<0)+180;    
end

[An1,v12,azm12] = AnisoFit(azmbin,va,azmo,360,1);
[An2,v22,azm22] = AnisoFit(azmbin,va,azmo,360,0);
[An3,v32,azm32] = AnisoFit(azmbin,va,azm2,180,0);


vfit22=fun1(An2.pout,(0:360)/180*pi);
vfit32=fun1(An3.pout,(0:360)/180*pi);

plotpar.st=st;
plotpar.evla=evla;
plotpar.evlo=evlo;
plotpar.azmo=azmo;
plotLoc(plotpar,handles,par) %-----plotLoc-----------


plotpar.azmf=azmf;
plotpar.vf=vf;
plotpar.azmfit=0:360;

plotpar.Aniso=An3;
plotpar.azma=azm2;
plotpar.va=va;
plotpar.azm2=azm32;
plotpar.v2=v32;
plotAnisoFit1(plotpar,handles.Plot1,par)%-----plotAnisoFit1-----------

plotpar.azm2=azm12;
plotpar.v2=v12;
plotpar.Aniso=An1;
plotpar.azma=azmo;
plotpar.va=va;
plotAnisoFit2(plotpar,handles.Plot2,par) %-----plotAnisoFit2-----------

plotpar.azm2=azm22;
plotpar.v2=v22;
plotpar.Aniso=An2;
plotpar.azma=azmo;
plotpar.va=va;
plotpar.vfit122=fun1(An1.pout(1:3),plotpar.azmfit/180*pi);
plotpar.vfit132=fun1(An3.pout,plotpar.azmfit/180*pi);
plotAnisoFit3(plotpar,handles.Plot3,par) %-----plotAnisoFit3-----------




function plotAnisoFit1(plotpar,Plotid,par)
Aniso=plotpar.Aniso;

azmf = plotpar.azmf;
vf = plotpar.vf;

azma = plotpar.azma;
va = plotpar.va;


azm2 = plotpar.azm2;
v2 = plotpar.v2;

azmfit=plotpar.azmfit;
vfit=fun1(Aniso.pout(1:3),azmfit/180*pi);
try
subplot(Plotid);
title(['(f) T=' num2str(par.period) ' s']) ;
catch
subplot(Plotid(1),Plotid(2),Plotid(3))
title(['(d) T=' num2str(par.period) ' s']) ;
end
plot(azmf,vf,'.','color',[0.5,0.5,0.5],'LineWidth',1.5);
hold on
scatter(azma,va,40,'.','MarkerEdgeColor',[0.1 0.5 0.8],'LineWidth',1.5);
scatter(azm2*180/pi,v2,60,'o','fill','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9 0.3 0.5])%
plot(azmfit,vfit,'color',[0.4660 0.7740 0.1880],'LineWidth',1.5);
N1=Aniso.viso-0.8;
k1=roundn(N1,-2) ;
N2=Aniso.viso+0.8;
k2=roundn(N2,-2) ;
set(gca,'ylim',[k1 k2]); %set(gca,'YTick',(k1:0.2:k2));
set(gca,'Xlim',[0 360]); set(gca,'XTick',0:60:360);
set(gca,'LineWidth',0.5,'FontSize',10);
text(10,k2-0.1,'mxazm:180, 2phase only')
text(10,k1+0.15,['MOA=' num2str(Aniso.M,'%.2f') '; FPD=' num2str(Aniso.Fai,'%.1f') '; v=' num2str(Aniso.viso,'%.2f')]) ;
title(['(d) T=' num2str(par.period) ' s']) ;
% xlabel('amz')
ylabel('v (km/s)') 
box on
hold off



function plotAnisoFit2(plotpar,Plotid,par)
st = plotpar.st;
evla = plotpar.evla;
evlo = plotpar.evlo;

azmf = plotpar.azmf;
vf = plotpar.vf;

azma = plotpar.azma;
va = plotpar.va;

Aniso = plotpar.Aniso;

v2 = plotpar.v2;
azm2 = plotpar.azm2;

azmfit=plotpar.azmfit;
vfit12=fun2(Aniso.pout,azmfit/180*pi);
vfit122=fun1(Aniso.pout(1:3),azmfit/180*pi);
vfit121=vfit12-vfit122+mean(vfit12);

try
subplot(Plotid);
title(['(e) T=' num2str(par.period) ' s']) ;
catch
subplot(Plotid(1),Plotid(2),Plotid(3))
title(['(d) T=' num2str(par.period) ' s']) ;
end

plot(azmf,vf,'.','color',[0,0,0],'LineWidth',1.5);
hold on
scatter(azma,va,5,'o','MarkerFaceColor',[0.1 0.5 0.8],'MarkerEdgeColor',[0.1 0.5 0.8],'LineWidth',1);
scatter(azm2*180/pi,v2,60,'o','fill','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9 0.3 0.5])%
h=plot((0:360),vfit121,'b-','LineWidth',1.5);
plot((0:360),vfit12,'color',[1 0 0],'LineStyle','-','LineWidth',1.5);
plot((0:360),vfit122,'color',[0.4660 0.7740 0.1880],'LineWidth',1.5);

N1=Aniso.viso-0.8;
k1=roundn(N1,-2) ;
N2=Aniso.viso+0.8;
k2=roundn(N2,-2) ;
set(gca,'ylim',[k1 k2]); %set(gca,'YTick',(k1:0.2:k2));
set(gca,'LineWidth',0.5,'FontSize',10);
set(gca,'Xlim',[0 360]); set(gca,'XTick',0:60:360);
text(10,k2-0.1,'mxazm:360, 1phase+2phase')
text(10,k1+0.15,['MOA=' num2str(Aniso.M,'%.2f') '; FPD=' num2str(Aniso.Fai,'%.1f') '; v=' num2str(Aniso.viso,'%.2f')]) ;
title(['(e) T=' num2str(par.period) ' s']) ;
ylabel('v (km/s)') 
box on
hold off


function plotAnisoFit3(plotpar,Plotid,par)
Aniso=plotpar.Aniso;

azmf = plotpar.azmf;
vf = plotpar.vf;

azma = plotpar.azma;
va = plotpar.va;


azm2 = plotpar.azm2;
v2 = plotpar.v2;

azmfit=plotpar.azmfit;
vfit=fun1(Aniso.pout(1:3),azmfit/180*pi);
try
    subplot(Plotid);
    plot(azmf,vf,'.','color',[0.5,0.5,0.5]*0,'LineWidth',1.5);
    hold on
    title(['(f) T=' num2str(par.period) ' s']) ;
catch
    subplot(Plotid(1),Plotid(2),Plotid(3))
    plot(azmf,vf,'.','color',[0.5,0.5,0.5]*0,'LineWidth',1.5);
    hold on
    title(['(d) T=' num2str(par.period) ' s']) ;
end

scatter(azma,va,40,'.','MarkerEdgeColor',[0.1 0.5 0.8],'LineWidth',1.5);
scatter(azm2*180/pi,v2,60,'o','fill','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9 0.3 0.5])%
plot(azmfit,vfit,'color',[0.4660 0.7740 0.1880],'LineWidth',1.5);
N1=Aniso.viso-0.8;
k1=roundn(N1,-2) ;
N2=Aniso.viso+0.8;
k2=roundn(N2,-2) ;
set(gca,'ylim',[k1 k2]); %set(gca,'YTick',(k1:0.2:k2));
set(gca,'Xlim',[0 360]); set(gca,'XTick',0:60:360);
set(gca,'LineWidth',0.5,'FontSize',10);
text(10,k2-0.1,'mxazm:360, 2phase only')
text(10,k1+0.15,['MOA=' num2str(Aniso.M,'%.2f') '; FPD=' num2str(Aniso.Fai,'%.1f') '; v=' num2str(Aniso.viso,'%.2f')]) ;
% xlabel('amz')
ylabel('v (km/s)') 
box on
hold off



function plotLoc(plotpar,handles,par)
st = plotpar.st;
evla = plotpar.evla;
evlo = plotpar.evlo;
azmo = plotpar.azmo;
try
subplot(handles.Plot4);
hold off
load coast
mx=long-st(1,2);
my=lat;
ind=mx>180|mx<-180|my>90|my<-90;
my(ind)=nan;
mx(ind)=nan;
plot(mx,my,'color',[0.5 0.5 0.5])
hold on
scatter(0,st(1),'^k','fill')
set(handles.Plot4,'ytick',[]);
set(handles.Plot4,'xtick',[]);
title('(a)')
plot(evlo-st(1,2),evla,'r.')
xlim([min(mx(~isnan(mx))) max(mx(~isnan(mx)))])
ylim([-90 90])
% ylim([max([min(my(~isnan(my))),-90]) min([90,max(my(~isnan(my)))])])
%%
subplot(handles.Plot6);
histogram(azmo,par.MaxAzm/par.azmbin)
% rose(azm0/180*pi)
xlim([-5 365])
xlabel('azm')
ylabel('num')
title('(b)')
box on

subplot(handles.Plot5);
title('(c)')
try
    delete(handles.Plot5.Children(1));
h=scatter (st(1,2),st(1,1),60,'Marker','o','MarkerEdgeColor',[1,1,1],'MarkerFaceColor',[0.2,0.2,0.7]);
    % axis equal
    xlim(par.Lonlim)
    ylim(par.Latlim)
catch
h=scatter (st(1,2),st(1,1),60,'Marker','o','MarkerEdgeColor',[1,1,1],'MarkerFaceColor',[0.2,0.2,0.7]);
    % axis equal
    xlim(par.Lonlim)
    ylim(par.Latlim)
end
box on


catch
subplot(2,2,handles(1));
hold off
mx=par.long-st(1,2);
my=par.lat;
ind=mx>180|mx<-180|my>90|my<-90;
my(ind)=nan;
mx(ind)=nan;
plot(mx,my,'color',[0.5 0.5 0.5])
hold on
scatter(0,st(1),'^k','fill')
set(gca,'ytick',[]);
set(gca,'xtick',[]);
xlabel('longitude')
ylabel('latitude')
title('(a)')
plot(evlo-st(1,2),evla,'r.')
xlim([min(mx(~isnan(mx))) max(mx(~isnan(mx)))])
ylim([-90 90])

% ylim([max([min(my(~isnan(my))),-90]) min([90,max(my(~isnan(my)))])])
%%
subplot(2,2,handles(2));
histogram(azmo,par.MaxAzm/par.azmbin)
% rose(azm0/180*pi)
xlim([-5 365])
xlabel('azm')
ylabel('num')
title('(b)')
box on

hs=subplot(2,2,handles(3));
title('(c)')
h=scatter (st(1,2),st(1,1),20,'Marker','.','MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[0.2,0.2,0.7]);
axis equal
xlim(par.Lonlim)
ylim(par.Latlim)
box on
end

% --- Executes during object creation, after setting all properties.
function stationList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stationList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PeriodList.
function PeriodList_Callback(hObject, eventdata, handles)
% hObject    handle to PeriodList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PeriodList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PeriodList
userdata=get(gcf,'userdata');             
par=userdata.par;
fout0= get(hObject,'String');
par.kinfile=fout0{get(hObject,'Value')};
par.period=str2double(par.kinfile);
try
    smr=get(handles.smr,'string');
    smr=strsplit(smr,',');
    for ism=1:length(smr)
        smris(ism)=str2double(cell2mat(smr(ism)));
    end
    if ~isempty(find(isnan(smris), 1))
        smr=get(handles.smr,'string');
        smr=strsplit(smr,' ');
        for ism=1:length(smr)
            smris(ism)=str2double(cell2mat(smr(ism)));
        end
    end
catch
end
pdb=[min(par.periods)-1,par.smrTbs,max(par.periods)+1];
smrs(length(par.periods))=0;
for i=1:length(pdb)-1
    id1= par.periods>pdb(i)&par.periods<=pdb(i+1)+1;
    smrs(id1)=smris(i);
end
par.smrs=smrs;
userdata.par=par;
set(gcf,'userdata',userdata);

% --- Executes during object creation, after setting all properties.
function PeriodList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PeriodList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AutoFitting.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------- AutoFitting_Callback ---------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AutoFitting_Callback(hObject, eventdata, handles)
% hObject    handle to AutoFitting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% get par
par=getparameters(handles);
is_parfor=par.is_parfor;
% get st list files
stas=dir([par.indir par.fs '*.mat']);
if isempty(stas)
    disp('no data')
end

par.stas=stas;
if isempty(par.outdir)||isempty(par.indir)
    disp('outdir ro indir is empty')
    return
end

% read faults
if par.is_figure==1 && is_parfor==0
    try
        par.faults=readfault(par.FF);
        par.Blocks=readfault(par.BBF);
        Blocks=par.Blocks;
        faults=par.faults;
        figure(1)          
        
        subplot(2,2,3)
        hold off
        for bi=1:length(Blocks)
            if ~isempty(Blocks(bi).a)
                plot(Blocks(bi).a(:,1), Blocks(bi).a(:,2),'color',[0.1660 0.540 0.1880],'Linewidth', 2)
                hold on
            end
        end
        
        for bi=1:length(faults)
            if ~isempty(faults(bi).a)
                plot(faults(bi).a(:,1), faults(bi).a(:,2),'k','Linewidth', 0.5)
            end
        end
        xlim(par.Lonlim)
        ylim(par.Latlim)
    catch
    end
end

para_aniso_parfor='para_aniso_parfor.mat';
save(para_aniso_parfor,'par');
disp('Aniso fitting...')
%% step 1: loops for grid points
 if is_parfor==1 
    parfor ist=1:length(stas)
        try
            parforFitAniso(ist,para_aniso_parfor)
        catch
            erf=fopen('errlog.txt','a+');
            fprintf(erf,'%f4.0\n',ist);
            fclose(erf);
            parforFitAniso(ist,para_aniso_parfor)
        end
    end
 elseif is_parfor==0  
     for ist=1:length(stas)
        parforFitAniso(ist,para_aniso_parfor)
    end
 end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------------------------------------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function parforFitAniso(kist,para_aniso_parfor)
load(para_aniso_parfor)
stas=par.stas; %#ok<*NODEF>
fs=par.fs;
indir= par.indir;
outdir = par.outdir;
is_figure=par.is_figure;
is_overwrite = par.is_overwrite;      
is_onephase = par.is_onephase;        
velocityType=par.vtype;   
is_parfor = par.is_parfor;         
is_bootstrap = [par.is_bootstrap par.bootTimes 0.9]; 
is_CIR=par.is_CIR;
periods = par.periods; 
vlim =par.vlim;    
azmban = par.azmbin;   
MaxAzm = par.MaxAzm;  
pbeg=min(par.periods)-1;
pend=max(par.periods)+1;
pdb=[pbeg par.smrTbs pend]; 
smrs=par.smrs;  
% smrb=par.smrs;  
% smrs(length(pdb)-1)=0;
% for i=1:length(pdb)-1
%     id1=par.periods>pdb(i)&par.periods<=pdb(i+1)+1;
%     smrs(id1)=smrb(i);
% end
%---------------------------------------------
bootNum = is_bootstrap(2);%
if ~exist(outdir,'dir') 
    mkdir(outdir)
end



temp = load([indir fs stas(kist).name]);
fnm=cell2mat(fieldnames(temp));
sti = getfield(temp,fnm); %#ok<GFLD>
Anisoi.st = [sti.st(1),sti.st(2)];

if strcmp(velocityType,'strc')
    V=[sti.vcr];
    Va=[sti.vcra];
    outfile = [outdir 'Aniso' num2str(MaxAzm) '_cr_' stas(kist).name];
elseif strcmp(velocityType,'dym')
    V=[sti.v];
    Va=[sti.va];
    outfile = [outdir 'Aniso' num2str(MaxAzm) '_' stas(kist).name];
end

if is_overwrite==0 && exist(outfile,'file')
    return
end

stPrd=[sti.periods];
st=sti.st;

evlas=sti.evla;
evlos=sti.evlo;
azmos=sti.azm;
dazms=sti.dazm;
disp(['ist = ' num2str(kist)])
for iprd = 1:length(periods)
    par.period=periods(iprd);
    id=find(stPrd==periods(iprd));
    if isempty(id)
        Anisoi.period(iprd)=periods(iprd);
        Anisoi.viso(iprd) = nan;
        Anisoi.stdv(iprd) = nan;
        Anisoi.Fai(iprd) = nan;
        Anisoi.M(iprd) = nan;
        Anisoi.a(iprd) = nan;
        Anisoi.b(iprd) = nan;
        Anisoi.fitRsrm(iprd) = nan;
        Anisoi.stdFai(iprd) = nan;
        Anisoi.stdM(iprd) = nan;
        Anisoi.stdA(iprd) = nan;
        Anisoi.stdB(iprd) = nan;
        continue
    end

    Anisoi.period(iprd)=periods(iprd);
    v = V(:,id);
    va =  Va(:,id);
    azma = azmos(:,id);
    dazm=dazms(:,id);
    evla=evlas;
    evlo=evlos;
    
    % remove bad velue
    ind=isnan(v) | isnan(va) |abs(dazm)>30 | v==0 | va==0;
    % ind=isnan(v);
    azma(ind)=[];
    v(ind)=[];
    va(ind)=[];
    dazm(ind)=[];
    evla(ind)=[];
    evlo(ind)=[];
    viso=median(va);
    
    % remove outlier 1
    ind=abs(va-viso)>vlim;
    azma(ind)=[];
    v(ind)=[];
    va(ind)=[];
    dazm(ind)=[];
    evla(ind)=[];
    evlo(ind)=[];
    viso=median(va);
    azma(azma<0)=azma(azma<0)+360;
    azma(azma>360)=azma(azma>360)-360;
    
    if length(v)<8&&~isempty(v)
        Anisoi.period(iprd)=periods(iprd);
        Anisoi.viso(iprd) = median(v);
        Anisoi.stdv(iprd) = nan;
        Anisoi.Fai(iprd) = nan;
        Anisoi.M(iprd) = nan;
        Anisoi.a(iprd) = nan;
        Anisoi.b(iprd) = nan;
        Anisoi.fitRsrm(iprd) = nan;
        Anisoi.stdFai(iprd) = nan;
        Anisoi.stdM(iprd) = nan;
        Anisoi.stdA(iprd) = nan;
        Anisoi.stdB(iprd) = nan;
        continue
    elseif isempty(v)
        Anisoi.period(iprd)=periods(iprd);
        Anisoi.viso(iprd) = nan;
        Anisoi.stdv(iprd) = nan;
        Anisoi.Fai(iprd) = nan;
        Anisoi.M(iprd) = nan;
        Anisoi.a(iprd) = nan;
        Anisoi.b(iprd) = nan;
        Anisoi.fitRsrm(iprd) = nan;
        Anisoi.stdFai(iprd) = nan;
        Anisoi.stdM(iprd) = nan;
        Anisoi.stdA(iprd) = nan;
        Anisoi.stdB(iprd) = nan;
        continue
    end
    
    % remove outlier 2
    if is_CIR==1
        An1 = AnisoFit(par.azmbin,va,azma,360,1);
        vfit=fun1(An1.pout,azma/180*pi);
        dv=v-vfit;
        ind=abs(dv)>vlim*0.5;
        azma(ind)=[];
        v(ind)=[];
        va(ind)=[];
        dazm(ind)=[];
        evla(ind)=[];
        evlo(ind)=[];
        viso=median(v);
    end
  
    if MaxAzm==180;
        while ~isempty(find(azma>180|azma<0,1))
            azma(azma>180) = azma(azma>180)-180;
            azma(azma<0) = azma(azma<0)+180;
        end
    end
    
    %%
    if length(v)<8&&~isempty(v)
        Anisoi.period(iprd)=periods(iprd);
        Anisoi.viso(iprd) = median(v);
        Anisoi.stdv(iprd) = nan;
        Anisoi.Fai(iprd) = nan;
        Anisoi.M(iprd) = nan;
        Anisoi.a(iprd) = nan;
        Anisoi.b(iprd) = nan;
        Anisoi.fitRsrm(iprd) = nan;
        Anisoi.stdFai(iprd) = nan;
        Anisoi.stdM(iprd) = nan;
        Anisoi.stdA(iprd) = nan;
        Anisoi.stdB(iprd) = nan;
        continue
    elseif isempty(v)
        Anisoi.period(iprd)=periods(iprd);
        Anisoi.viso(iprd) = nan;
        Anisoi.stdv(iprd) = nan;
        Anisoi.Fai(iprd) = nan;
        Anisoi.M(iprd) = nan;
        Anisoi.a(iprd) = nan;
        Anisoi.b(iprd) = nan;
        Anisoi.fitRsrm(iprd) = nan;
        Anisoi.stdFai(iprd) = nan;
        Anisoi.stdM(iprd) = nan;
        Anisoi.stdA(iprd) = nan;
        Anisoi.stdB(iprd) = nan;
        continue
    end

    %%
    bandNum = MaxAzm/azmban;    
    ig=0;
    for ibn = 1:bandNum
        bd1 = (ibn-1)*azmban; bd2 = ibn*azmban;
        vID = find(azma>=bd1 & azma<=bd2, 1);
        if isempty(vID)
            continue
        else
            ig=ig+1;
        end
    end
    %%
    if  ig<4
        Anisoi.viso(iprd) = median(v);
        Anisoi.stdv(iprd) = nan;
        Anisoi.Fai(iprd) = nan;
        Anisoi.M(iprd) = nan;
        Anisoi.a(iprd) = nan;
        Anisoi.b(iprd) = nan;
        Anisoi.fitRsrm(iprd) = nan;
        Anisoi.stdFai(iprd) = nan;
        Anisoi.stdM(iprd) = nan;
        Anisoi.stdA(iprd) = nan;
        Anisoi.stdB(iprd) = nan;
        continue
    end
    %% (azmban,v,azm,MaxAzm,is_onephase)
    [An,v2,azm2] = AnisoFit(azmban,va,azma,MaxAzm,is_onephase);
    Anisoi.viso(iprd) = median(v);
    Anisoi.Fai(iprd) = An.Fai;
    Anisoi.M(iprd) = An.M;
    Anisoi.fitRsrm(iprd) = An.fitRsrm;
    Anisoi.a(iprd) = An.a;
    Anisoi.b(iprd) = An.b;
    if is_bootstrap(1)==1
        fitNum = round(length(v)*is_bootstrap(3));
        bootAniso = bootaniso(fitNum,v,va,azma,azmban,bootNum,MaxAzm,is_onephase,vlim);
        bootFai=[bootAniso.Fai];        
        bootM=[bootAniso(:).Fai];
        boota=[bootAniso(:).a];
        bootb=[bootAniso(:).b];        
        nanID=find(isnan(bootFai));
        if length(nanID)/length(bootFai)>0.5
            Anisoi.stdFai(iprd) = nan;
            Anisoi.stdM(iprd) = nan;
            Anisoi.stdA(iprd) = nan;
            Anisoi.stdB(iprd) = nan;
            Anisoi.stdv(iprd) = nan;
            continue
        end
        bootFai(nanID)=[];
        bootM(nanID)=[];
        boota(nanID)=[];
        bootb(nanID)=[];
        dFai=bootFai-An.Fai;
        while ~isempty(find(abs(dFai)>90,1))
            bootFai(abs(dFai)>90)=bootFai(abs(dFai)>90)-180.*dFai(abs(dFai)>90)./abs(dFai(abs(dFai)>90));
            dFai=bootFai-median(bootFai);
        end
        Anisoi.stdv(iprd)=std([bootAniso.viso]);
        Anisoi.stdFai(iprd) = std(bootFai);
        Anisoi.stdM(iprd) = std([bootAniso.M]);
        Anisoi.stdA(iprd) = std([bootAniso.a]);
        Anisoi.stdB(iprd) = std([bootAniso.b]);
    else
        
        Anisoi.stdv(iprd)=nan;
        Anisoi.stdFai(iprd) = nan;
        Anisoi.stdM(iprd) = nan;
        Anisoi.stdA(iprd) = nan;
        Anisoi.stdB(iprd) = nan;
    end
    
    if par.is_figure==1 && par.is_parfor==0
        figure(1)
        plotpar.st=st;
        plotpar.evla=evla;
        plotpar.evlo=evlo;
        plotpar.azmo=azma;        
        plotLoc(plotpar,[1,2,3],par) %-----plotLoc-----------
        plotpar.azmf=azmos(:,id);
        plotpar.vf=Va(:,id);
        plotpar.azmfit=0:360;
        
        plotpar.Aniso=An;
        plotpar.azma=azma;
        plotpar.va=va;
        plotpar.azm2=azm2;
        plotpar.v2=v2;
        figure(1)
        if par.MaxAzm==180
            plotAnisoFit1(plotpar,[2,2,4],par)%-----plotAnisoFit1-----------
        elseif par.MaxAzm==360 && par.is_onephase==1
            plotAnisoFit2(plotpar,[2,2,4],par)%-----plotAnisoFit2-----------
        elseif par.MaxAzm==360 && par.is_onephase==0
            plotAnisoFit3(plotpar,[2,2,4],par)%-----plotAnisoFit3-----------
        end
        drawnow
    end
    
end
save(outfile,'Anisoi')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------- bootstrap estimate ---------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Aniso = bootaniso(fitNum,v,va,azmo,azmban,bootNum,MaxAzm,is_onephase,vlim)
%bootstrap method
bootID = 1:fitNum;
bootID = bootstrp(bootNum, @bootline,bootID);
% [Aniso(1:bootNum).flag]=deal(1);
for iboot = 1:bootNum
    bootVa=[];
    bootV=[];
    bootazm=[];
    ID=unique(bootID(iboot,:));
    bootVa(:,1) = va(ID);
    bootV(:,1) = v(ID);
    bootazm(:,1) = azmo(ID);
    bootdat = [bootVa,bootazm];
    bootdat=sortrows(bootdat,2);
    bootVa = bootdat(:,1);
    bootazm = bootdat(:,2);
    bandNum = MaxAzm/azmban;
    vt2(1:bandNum) = nan;
    for ibn = 1:bandNum
        bd1 = (ibn-1)*azmban; bd2 = ibn*azmban;
        vID = find(bootazm>=bd1 & bootazm<=bd2);
        if isempty(vID)
            continue
        else
            vt2(ibn) = median(bootVa(vID));
        end
    end
    vt=vt2(~isnan(vt2));
    if length(vt)<5
         Aniso(iboot).viso=median(bootV);
         Aniso(iboot).Fai=nan;
         Aniso(iboot).M=nan;
         Aniso(iboot).a=nan;
         Aniso(iboot).b=nan;
         Aniso(iboot).Fai1=nan;
         Aniso(iboot).M1=nan;
         Aniso(iboot).a1=nan;
         Aniso(iboot).b1=nan;
         Aniso(iboot).fitRsrm=nan;
         Aniso(iboot).pout=[nan,nan,nan];
         continue
    end
    Aniso(iboot) = AnisoFit(azmban,bootVa,bootazm,MaxAzm,is_onephase);
    Aniso(iboot).viso=median(bootV);
end




function idx=nrGrd(st,is,dsmax,xs,ys)
yi=st(is,1);
xi=st(is,2);
dx=abs(xs-xi);
dy=abs(ys-yi);
stid1=find(dx<4*dsmax & dy<4*dsmax);
xn=xs(stid1);
yn=ys(stid1);

dst=ones(size(xn))*nan;
for ist=1:length(xn)
    dst(ist)=distance(yi,xi,yn(ist),xn(ist));
end

idx=dst<=dsmax;
idx=stid1(idx);








%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------------------- Anisotropy fitting ---------------------
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Aniso,v2,azm2] = AnisoFit(azmban,v,azm,MaxAzm,is_onephase)
%          azmbin averege
bandNum = MaxAzm/azmban;
k=0;
for ibn = 1:bandNum
    bd1 = (ibn-1)*azmban;
    bd2 = ibn*azmban;
    vID = find(azm>=bd1 & azm<=bd2);
    if length(vID)<1
        continue
    else
        k = k+1;
        v2(k) = median(v(vID)); %#ok<*AGROW>
        azmb = azm(vID);
        ddv = abs(v(vID)-v2(k));
        azmID = ddv==min(ddv);
        azm2 (k) = mean(azmb(azmID));    
    end
end
viso = median(v);
azm2 = azm2/180*pi;
if k<5
    Aniso.viso = viso;
    Aniso.Fai = 0;
    Aniso.M = 0;
    Aniso.a = 0;
    Aniso.b = 0;
    Aniso.fitRsrm = nan;
    Aniso.Fai1 = 0;
    Aniso.M1 = 0;
    Aniso.a1 = 0;
    Aniso.b1 = 0; 
    Aniso.pout=[viso,0,0,0,0];
    disp('no enouth data, skip');
    return
end

if is_onephase==0 || MaxAzm==180
    lb=[viso*0.8,-0.15,-0.15];
    ub=[viso*1.2,0.15,0.15];
    co = [viso,0,0];
    options = optimoptions('lsqcurvefit','Display','off');
    [pout,anisoRsrm] = lsqcurvefit(@fun1,co,azm2,v2,lb,ub,options);
    x=0:MaxAzm;
    y=fun1(pout,x/180*pi);
    y1=y*0+viso;
    FPD=x(y==max(y));
    FPD=FPD(1);
    MOA=(max(y)-min(y))/viso*100;
    FPD(FPD>180)=FPD(FPD>180)-180;
    Aniso.viso = pout(1);
    Aniso.Fai = FPD;
    Aniso.M = MOA;
    Aniso.a = pout(2);
    Aniso.b = pout(3);
    Aniso.fitRsrm = anisoRsrm;
    Aniso.Fai1 = nan;
    Aniso.M1 = nan;
    Aniso.a1 = nan;
    Aniso.b1 = nan;
else
    lb=[viso*0.8,-0.15,-0.15,-0.15,-0.15];
    ub=[viso*1.2,0.15,0.15,0.15,0.15];
    co = [viso,0,0,0,0];
    options = optimoptions('lsqcurvefit','Display','off');
    [pout,anisoRsrm] = lsqcurvefit(@fun2,co,azm2,v2,lb,ub,options);
    x=0:360;
    y=fun1(pout(1:3),x/180*pi);
    y1=fun2(pout,x/180*pi)-y+viso;
    FPD=x(y==max(y));
    FPD=FPD(1);
    FPD(FPD>180)=FPD(FPD>180)-180;
    MOA=(max(y)-min(y))/viso*100;
    FPD(FPD>180)=FPD(FPD>180)-180;
    MOA1=(max(y1)-min(y1))/viso*100;
    FPD1=x(y1==max(y1));
    Aniso.viso = pout(1);
    Aniso.Fai = FPD;
    Aniso.M = MOA;
    Aniso.a = pout(2);
    Aniso.b = pout(3);
    Aniso.fitRsrm = anisoRsrm;
    Aniso.Fai1 = FPD1;
    Aniso.M1 = MOA1;
    Aniso.a1 = pout(4);
    Aniso.b1 = pout(5);
end
Aniso.pout=pout;


function y1=bootline(y)
y1=y;

function fy1=fun1(abc,fx)
fy1=abc(1)+abc(2)*cos(2*fx)+abc(3)*sin(2*fx);


function fy1=fun2(abc,fx)
fy1=abc(1)+abc(2)*cos(2*fx)+abc(3)*sin(2*fx)+abc(4)*cos(fx)+abc(5)*sin(fx);





function indir_Callback(hObject, eventdata, handles)
% hObject    handle to indir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of indir as text
%        str2double(get(hObject,'String')) returns contents of indir as a double
userdata=get(gcf,'userdata');             
par=userdata.par;
par.indir=get(handles.indir,'string');
userdata.par=par;
set(gcf,'userdata',userdata);


% --- Executes during object creation, after setting all properties.
function indir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to indir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SLTindir.
function SLTindir_Callback(hObject, eventdata, handles)
% hObject    handle to SLTindir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');             
par=userdata.par;
dname = uigetdir('./');
par.indir=[dname par.fs];

if dname==0
    return
end

set(handles.indir,'string',par.indir)
stas=dir([par.indir par.fs '*.mat']);
if isempty(stas)
    disp('no data')
end
par.stas=stas;
pstr(1:length(par.periods))={nan};
set(handles.stationList,'string',{stas.name});

userdata.par=par;
set(gcf,'userdata',userdata);


function outdir_Callback(hObject, eventdata, handles)
% hObject    handle to outdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outdir as text
%        str2double(get(hObject,'String')) returns contents of outdir as a double
userdata=get(gcf,'userdata');             
par=userdata.par;
par.outdir=get(handles.outdir,'string');
userdata.par=par;
set(gcf,'userdata',userdata);


% --- Executes during object creation, after setting all properties.
function outdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SLToutdir.
function SLToutdir_Callback(hObject, eventdata, handles)
% hObject    handle to SLToutdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');             
par=userdata.par;
dname = uigetdir('./');
if dname==0
    return
end
par.outdir=[dname par.fs];
set(handles.outdir,'string',par.outdir)
userdata.par=par;
set(gcf,'userdata',userdata);


% --- Executes on button press in MergeAniso.
function MergeAniso_Callback(hObject, eventdata, handles)
% hObject    handle to MergeAniso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');             
par=userdata.par;
if isempty(par.indir) || isempty(par.outdir)
    disp('outdir or outdir2 is empty')
    return
end
MaxAzm=par.MaxAzm;
indir=par.outdir;
vtype=par.vtype;
fs=par.fs;

if strcmp(vtype,'strc')
    anisoDirs = dir([indir fs 'Aniso' num2str(MaxAzm) '_cr_bk_*_*.mat']);
    outmatf=[par.outdir2 fs 'Aniso' num2str(MaxAzm) 'azmbin' num2str(par.azmbin) 'onephs' num2str(par.is_onephase) '_cr_allst.mat'];
elseif strcmp(vtype,'dym')
    anisoDirs = dir([indir fs 'Aniso' num2str(MaxAzm) '_bk_*_*.mat']);
    outmatf=[par.outdir2 fs 'Aniso' num2str(MaxAzm) 'azmbin' num2str(par.azmbin) 'onephs' num2str(par.is_onephase) '_allst.mat'];
end
par.outmatf=outmatf;
% par.outmatfs=outmatfs;
if length(anisoDirs)>1
    disp('mergering aniso to a noe file... ')
    disp(outmatf)
else
    disp('no data');
    return
end

for ia =1:length(anisoDirs)
    temp = load([indir fs anisoDirs(ia).name]);
    fnm=cell2mat(fieldnames(temp));
    Anisoi = getfield(temp,fnm); %#ok<GFLD>
    Aniso.st(ia,:)=Anisoi.st;
    Aniso.viso(ia,:) = Anisoi.viso;
    Aniso.stdv(ia,:) = Anisoi.stdv;
    Aniso.Fai(ia,:) =  Anisoi.Fai;
    Aniso.M(ia,:) = Anisoi.M;
    Aniso.a(ia,:) = Anisoi.a;
    Aniso.b(ia,:) = Anisoi.b;
    Aniso.fitRsrm(ia,:) = Anisoi.fitRsrm;
    Aniso.stdFai(ia,:) = Anisoi.stdFai;
    Aniso.stdM(ia,:) = Anisoi.stdM;
    Aniso.stdA(ia,:) = Anisoi.stdA;
    Aniso.stdB(ia,:) = Anisoi.stdB;

    stla=Anisoi.st(1);
    stlo=Anisoi.st(2);
    v=Anisoi.viso;
    dv=Anisoi.stdv;
    FPD=Anisoi.Fai;
    dfpd=Anisoi.stdFai;
    MOA=Anisoi.M;
    dmoa=Anisoi.stdM;
    a=Anisoi.a;
    b=Anisoi.b;
    period=Anisoi.period;
    nanID=isnan(MOA);
    
    v(nanID)=[];
    dv(nanID)=[];
    FPD(nanID)=[];
    dfpd(nanID)=[];
    MOA(nanID)=[];
    dmoa(nanID)=[];
    a(nanID)=[];
    b(nanID)=[];
    period(nanID)=[];
%     if ~isempty(period)
%         for j=1:length(period)
%             pdi=period(j);
%             outtxtf=[outpath 'Aniso' num2str(MaxAzm) 'p' num2str(pdi) 'allst.txt'];
%             ouf=fopen(outtxtf,'a+');                
%             fprintf(ouf,'%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %5.0f\n',...
%                 stla,stlo,v(j),dv(j),FPD(j),dfpd(j),MOA(j),dmoa(j),a(j),b(j),pdi);
%             fclose(ouf);
%         end
%     end
end
Aniso.period = Anisoi.period;
save(outmatf,'Aniso');
disp(['save to: ' outmatf])
userdata.par=par;
set(gcf,'userdata',userdata);

% --- Executes on button press in SMOOTH.
function SMOOTH_Callback(hObject, eventdata, handles)
% hObject    handle to SMOOTH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
par=getparameters(handles);
[infile,path,indx] = uigetfile('*.mat');
strf=strsplit(infile,'.');
strf=cell2mat(strf(1));
if indx==0;
    return
end
outmatf=[path infile];
outmatf2=[par.outdir2 par.fs strf '_smr.mat'];
outpath=par.outdir2;
smrs=par.smrs;
MaxAzm=par.MaxAzm;
periods=par.periods;
is_overwrite=par.is_overwrite;
is_parfor=par.is_parfor;
para_aniso_parfor='para_aniso_parfor.mat';
load(outmatf)
Anisoin=Aniso; 
save(para_aniso_parfor,'smrs','outpath','MaxAzm','periods','Anisoin','is_overwrite');
if is_parfor==1
    parfor ipd=1:length(periods)        
        smoothaniso(ipd,para_aniso_parfor);
    end
elseif is_parfor==0
    for ipd=1:length(periods)        
        smoothaniso(ipd,para_aniso_parfor);
    end
end
outmats=dir([outpath par.fs 'Aniso' num2str(MaxAzm) 'p*allst_smr.mat']);
for ipd=1:length(outmats)
    Anisoi=load([outpath par.fs outmats(ipd).name]);
    Anisoi=Anisoi.Aniso;
    pdx=periods==Anisoi.period;
    Aniso.viso(:,pdx) = Anisoi.viso;
    Aniso.Fai(:,pdx) = Anisoi.Fai;
    Aniso.M(:,pdx) = Anisoi.M;
    Aniso.a(:,pdx) = Anisoi.a;
    Aniso.b(:,pdx) = Anisoi.b;
    Aniso.stdv(:,pdx) = Anisoi.stdv;
    Aniso.stdFai(:,pdx) = Anisoi.stdFai;
    Aniso.stdM(:,pdx) = Anisoi.stdM;
    Aniso.stdA(:,pdx) = Anisoi.stdA;
    Aniso.stdB(:,pdx) = Anisoi.stdB;
end
save(outmatf2,'Aniso')


function smoothaniso(ipd,parmat)
disp(ipd)
load(parmat)
pdi=periods(ipd);
smr=smrs(periods==pdi);
if isempty(smr)||smr==0
    return
end
st=Anisoin.st;
[nst,~]=size(st);
xo=st(:,2);
ya=st(:,1);
outmat=[outpath 'Aniso' num2str(MaxAzm) 'p' num2str(pdi) 'allst_smr.mat'];
outtxf=[outpath 'Aniso' num2str(MaxAzm) 'p' num2str(pdi) 'allst_smr.txt'];

if is_overwrite==0 && exist(outmat,'file')
    return
end

Vip=Anisoin.viso(:,ipd);
Aip=Anisoin.a(:,ipd);
Bip=Anisoin.b(:,ipd);
MOAip=Anisoin.M(:,ipd);
FPDip=Anisoin.Fai(:,ipd);
dvip=Anisoin.stdv(:,ipd);
dFPDip=Anisoin.stdFai(:,ipd);
dMOAip=Anisoin.stdM(:,ipd);
dAip=Anisoin.stdA(:,ipd);
dBip=Anisoin.stdB(:,ipd);
Aniso.st=st;
Aniso.period=pdi;
ouf=fopen(outtxf,'w');
for ist=1:nst
    idx=nrGrd(st,ist,smr,xo,ya);
    idx(isnan(Aip(idx))|dFPDip(idx)>45|dMOAip(idx)>median(MOAip(~isnan(MOAip)|MOAip>0))|MOAip(idx)<=0)=[];
    if length(idx)<1
        Aniso.viso(ist,1)=nan;
        Aniso.Fai(ist,1) = 0;
        Aniso.M(ist,1) = 0;
        Aniso.a(ist,1) = 0;
        Aniso.b(ist,1) = 0;
        Aniso.stdv(ist,1)=nan;
        Aniso.stdFai(ist,1) = nan;
        Aniso.stdM(ist,1) = nan;
        Aniso.stdA(ist,1) = nan;
        Aniso.stdB(ist,1) = nan;
        continue
    end
    MOAis=MOAip(idx);
    difM=MOAis/mean(MOAis);
    idx(difM>1.8|difM<0.2)=[];
    if length(idx)<1
        Aniso.viso(ist,1)=nan;
        Aniso.Fai(ist,1) = 0;
        Aniso.M(ist,1) = 0;
        Aniso.a(ist,1) = 0;
        Aniso.b(ist,1) = 0;
        Aniso.stdv(ist,1)=nan;
        Aniso.stdFai(ist,1) = nan;
        Aniso.stdM(ist,1) = nan;
        Aniso.stdA(ist,1) = nan;
        Aniso.stdB(ist,1) = nan;
        continue
    end
    vis=Vip(idx);
    Ais=Aip(idx);
    Bis=Bip(idx);
    dvis=dvip(idx);
    dAis=dAip(idx);
    dBis=dBip(idx);
    dFPDis=dFPDip(idx);
    dMOAis=dMOAip(idx);    
    azm=(0:360);
    vazm=fun1([mean(vis),mean(Ais),mean(Bis)],azm*pi/180);
    plot(vazm,'r');
    hold on

    cofR=zeros(size(vis));
    for i=1:length(vis)
        vazmi=fun1([vis(i),Ais(i),Bis(i)],azm*pi/180);
        plot(vazmi,'k');
        cofRi=corrcoef(vazmi,vazm);
        cofR(i)=cofRi(1,2);
    end
    idx2=cofR>=0.5;
    if isempty(find(idx2,1))
        Aniso.viso(ist,1)=mean(Vip(idx));
        Aniso.Fai(ist,1) = 0;
        Aniso.M(ist,1) = 0;
        Aniso.a(ist,1) = 0;
        Aniso.b(ist,1) = 0;
        Aniso.stdv(ist,1)=nan;
        Aniso.stdFai(ist,1) = nan;
        Aniso.stdM(ist,1) = nan;
        Aniso.stdA(ist,1) = nan;
        Aniso.stdB(ist,1) = nan;
        continue
    end
    
    
    
    vis_m=mean(vis(idx2,:));
    Ais_m=mean(Ais(idx2,:));
    Bis_m=mean(Bis(idx2,:));
    
    dAis_m=mean(dAis(idx2,:));    
    dv_m=mean(dvis(idx2,:));
    dBis_m=mean(dBis(idx2,:));
    dFPDis_m=mean(dFPDis(idx2,:));
    dMOAis_m=mean(dMOAis(idx2,:));
    
    vazm=fun1([mean(vis(idx2,:)),Ais_m,Bis_m],azm*pi/180);
    FPD=azm(vazm==max(vazm));
    FPD=FPD(1);
    FPD(FPD>180)=FPD(FPD>180)-180;
    MOA=(max(vazm)-min(vazm))/mean(vazm)*100;
    FPD(FPD>180)=FPD(FPD>180)-180;
    
    Aniso.viso(ist,1) = vis_m;
    Aniso.Fai(ist,1) = FPD;
    Aniso.M(ist,1) = MOA;
    Aniso.a(ist,1) = Ais_m;
    Aniso.b(ist,1) = Bis_m;
    Aniso.stdv(ist,1) = dv_m;
    Aniso.stdFai(ist,1) = dFPDis_m;
    Aniso.stdM(ist,1) = dMOAis_m;
    Aniso.stdA(ist,1) = dAis_m;
    Aniso.stdB(ist,1) = dBis_m;
    stla=ya(ist,1);
    stlo=xo(ist,1);
    if ~isnan(vis_m)
        fprintf(ouf,'%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %5.0f\n',...
            stla,stlo,vis_m,dv_m,FPD,dFPDis_m,MOA,dMOAis_m,Ais_m,Bis_m,pdi);
    end
    
end
fclose(ouf);
save(outmat,'Aniso')


function outdir2_Callback(hObject, eventdata, handles)
% hObject    handle to outdir2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outdir2 as text
%        str2double(get(hObject,'String')) returns contents of outdir2 as a double
userdata=get(gcf,'userdata');
par=userdata.par;
par.outdir2=get(handles.outdir2,'string');
userdata.par=par;
set(gcf,'userdata',userdata);

% --- Executes during object creation, after setting all properties.
function outdir2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outdir2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SLToutdir2.
function SLToutdir2_Callback(hObject, eventdata, handles)
% hObject    handle to SLToutdir2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');             
par=userdata.par;
dname = uigetdir('./');
par.outdir2=[dname par.fs];
set(handles.outdir2,'string',par.outdir2)
userdata.par=par;
set(gcf,'userdata',userdata);

% --- Executes on button press in PLOTMAP.
function PLOTMAP_Callback(hObject, eventdata, handles) %#ok<*INUSL>
% hObject    handle to PLOTMAP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
par=getparameters(handles);
try
    [file,path,indx] = uigetfile([par.outdir2 '*.mat']);
catch
    [file,path,indx] = uigetfile('*.mat');
end
if indx==0;
    return
end
load([path file]);
x=Aniso.st(:,2);
y=Aniso.st(:,1);
xg=unique(x);
yg=unique(y);
[xg,yg] = meshgrid(xg,yg);
try
    for i=par.period
        periods=unique(Aniso.period);
        ids=find(periods==i);
        Visos=Aniso.viso(:,ids);
        FPD=Aniso.Fai(:,ids);
        MOA=Aniso.M(:,ids)/100;
        MOA(MOA>5)=5;
        
        v1=median(Visos(~isnan(Visos)))-4*std(Visos(~isnan(Visos)));
        v2=median(Visos(~isnan(Visos)))+4*std(Visos(~isnan(Visos)));
        deleID=Visos>v2|Visos<v1;
        
        
        Visos(deleID)=nan;
        F(deleID)=nan;
        Visos(deleID)=nan;
        
        a=MOA.*sin(FPD/180*pi);
        b=MOA.*cos(FPD/180*pi);
        [n,m]=size(xg);
        for isg=1:n
            for jsg=1:m
                idx=find(x==xg(isg,jsg) & y==yg(isg,jsg));
                if isempty(idx)
                    vGrd(isg,jsg)=nan;
                    ag(isg,jsg)=nan;
                    bg(isg,jsg)=nan;
                else
                    vGrd(isg,jsg)=Visos(idx);
                    ag(isg,jsg)=a(idx);
                    bg(isg,jsg)=b(idx);
                end
                
                
            end
        end
        figure
%         subplot(handles.Plot5);
        vGrd=griddata(x,y,Visos,xg,yg);
        hold off
        [~,c]=contourf(xg,yg,vGrd,200);
        colormap(flipud(jet))
        set(c,'edgecolor','none');
        % c.LineColor='w';
        hold on        
        
        faults=readfault(par.FF);
        Blocks=readfault(par.BBF);  
        for bi=1:length(Blocks)
            if ~isempty(Blocks(bi).a)
                plot(Blocks(bi).a(:,1), Blocks(bi).a(:,2),'color',[0.1660 0.540 0.1880],'Linewidth', 2)
                hold on
            end
        end
        
        for bi=1:length(faults)
            if ~isempty(faults(bi).a)
                plot(faults(bi).a(:,1), faults(bi).a(:,2),'w','Linewidth', 0.5)
            end
        end
        
        
        h1=quiver(xg,yg,ag,bg,'k');
        h1.ShowArrowHead = 'off';
        h2=quiver(xg,yg,-ag,-bg,'k');
        h2.ShowArrowHead = 'off';
        axis equal
        title(['T = ' num2str(i) 's'])
        xlim(par.Lonlim)
        ylim(par.Latlim)
        hold off
    end
catch
    for i=Aniso.period
        periods=unique(Aniso.period);
        ids=find(periods==i);
        Visos=Aniso.viso(:,ids);
        FPD=Aniso.Fai(:,ids);
        MOA=Aniso.M(:,ids)/100;
        MOA(MOA>5)=5;
        
        v1=median(Visos(~isnan(Visos)))-4*std(Visos(~isnan(Visos)));
        v2=median(Visos(~isnan(Visos)))+4*std(Visos(~isnan(Visos)));
        deleID=Visos>v2|Visos<v1;
        
        
        Visos(deleID)=nan;
        F(deleID)=nan;
        Visos(deleID)=nan;
        
        a=MOA.*sin(FPD/180*pi);
        b=MOA.*cos(FPD/180*pi);
        [n,m]=size(xg);
        for isg=1:n
            for jsg=1:m
                idx=find(x==xg(isg,jsg) & y==yg(isg,jsg));
                if isempty(idx)
                    vGrd(isg,jsg)=nan;
                    ag(isg,jsg)=nan;
                    bg(isg,jsg)=nan;
                else
                    vGrd(isg,jsg)=Visos(idx);
                    ag(isg,jsg)=a(idx);
                    bg(isg,jsg)=b(idx);
                end
                
                
            end
        end
        figure
%         subplot(handles.Plot5);
        vGrd=griddata(x,y,Visos,xg,yg);
        hold off
        [~,c]=contourf(xg,yg,vGrd,200);
        colormap(flipud(jet))
        set(c,'edgecolor','none');
        hold on
        
        faults=readfault(par.FF);
        Blocks=readfault(par.BBF);
        for bi=1:length(Blocks)
            if ~isempty(Blocks(bi).a)
                plot(Blocks(bi).a(:,1), Blocks(bi).a(:,2),'color',[0.1660 0.540 0.1880],'Linewidth', 2)
                hold on
            end
        end
        
        for bi=1:length(faults)
            if ~isempty(faults(bi).a)
                plot(faults(bi).a(:,1), faults(bi).a(:,2),'w','Linewidth', 0.5)
            end
        end
        
        h1=quiver(xg,yg,ag,bg,'k');
        h1.ShowArrowHead = 'off';
        h2=quiver(xg,yg,-ag,-bg,'k');
        h2.ShowArrowHead = 'off';
        axis equal
        title(['T = ' num2str(i) 's'])
        xlim(par.Lonlim)
        ylim(par.Latlim)
        hold off
    end
end

% --- Executes on button press in is_overwrite.
function is_overwrite_Callback(hObject, eventdata, handles)
% hObject    handle to is_overwrite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of is_overwrite
userdata=get(gcf,'userdata');
par=userdata.par;
par.is_overwrite=get(handles.is_overwrite,'value');
userdata.par=par;
set(gcf,'userdata',userdata);

% --- Executes on button press in faultf.
function faultf_Callback(hObject, eventdata, handles)
% hObject    handle to faultf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par=userdata.par;
[file,path,indx] = uigetfile('*.txt');
if indx==0;
    return
end
par.FF=[path file];
try
faults=readfault(par.FF);
par.faults=faults;

subplot(handles.Plot5);
for bi=1:length(faults)
    if ~isempty(faults(bi).a)
        plot(faults(bi).a(:,1), faults(bi).a(:,2),'k','Linewidth', 0.5)
        hold on
    end
end
xlim(par.Lonlim)
ylim(par.Latlim)
% axis equal
catch
end
userdata.par=par;
set(gcf,'userdata',userdata);

% --- Executes on button press in BBF.
function BBF_Callback(hObject, eventdata, handles)
% hObject    handle to BBF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par=userdata.par;
[file,path,indx] = uigetfile('*.txt');
if indx==0;
    return
end
par.BBF=[path file];

try
Blocks=readfault(par.BBF);
par.Blocks=Blocks;
subplot(handles.Plot5);
for bi=1:length(Blocks)
    if ~isempty(Blocks(bi).a)
        plot(Blocks(bi).a(:,1), Blocks(bi).a(:,2),'color',[0.1660 0.540 0.1880],'Linewidth', 2)
        hold on
    end
end
xlim(par.Lonlim)
ylim(par.Latlim)
% axis equal
catch
end

userdata.par=par;
set(gcf,'userdata',userdata);


function xylim_Callback(hObject, eventdata, handles)
% hObject    handle to xylim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xylim as text
%        str2double(get(hObject,'String')) returns contents of xylim as a double
userdata=get(gcf,'userdata');
par=userdata.par;
xylim=get(handles.xylim,'string');
try
    xylim=strsplit(xylim,',');
    Latlim=[str2double(cell2mat(xylim(1))),str2double(cell2mat(xylim(2)))];
    Lonlim=[str2double(cell2mat(xylim(3))),str2double(cell2mat(xylim(4)))];
catch
    xylim=strsplit(xylim,' ');
    Latlim=[str2double(cell2mat(xylim(1))),str2double(cell2mat(xylim(2)))];
    Lonlim=[str2double(cell2mat(xylim(3))),str2double(cell2mat(xylim(4)))];
end
par.Latlim=Latlim;
par.Lonlim=Lonlim;
userdata.par=par;
set(gcf,'userdata',userdata);

% --- Executes during object creation, after setting all properties.
function xylim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xylim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function smrTb_Callback(hObject, eventdata, handles)
% hObject    handle to smrTb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of smrTb as text
%        str2double(get(hObject,'String')) returns contents of smrTb as a double
userdata=get(gcf,'userdata');             
par=userdata.par;
try
    smrTb=get(handles.smrTb,'string');
    smrTb=strsplit(smrTb,',');
    for ism=1:length(smrTb)
        smrTbs(ism)=str2double(cell2mat(smrTb(ism)));
    end
    if ~isempty(find(isnan(smrTbs), 1))
        smrTb=get(handles.smrTb,'string');
        smrTb=strsplit(smrTb,' ');
        for ism=1:length(smrTb)
            smrTbs(ism)=str2double(cell2mat(smrTb(ism)));
        end
    end
catch
end
par.smrTbs=smrTbs;
userdata.par=par;
set(gcf,'userdata',userdata);

% --- Executes during object creation, after setting all properties.
function smrTb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smrTb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function smr_Callback(hObject, eventdata, handles)
% hObject    handle to smr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of smr as text
%        str2double(get(hObject,'String')) returns contents of smr as a double
userdata=get(gcf,'userdata');
par=userdata.par;

try
    smr=get(handles.smr,'string');
    smr=strsplit(smr,',');
    for ism=1:length(smr)
        smris(ism)=str2double(cell2mat(smr(ism)));
    end
    if ~isempty(find(isnan(smris), 1))
        smr=get(handles.smr,'string');
        smr=strsplit(smr,' ');
        for ism=1:length(smr)
            smris(ism)=str2double(cell2mat(smr(ism)));
        end
    end
catch
end
pdb=[min(par.periods)-1,par.smrTbs,max(par.periods)+1];
smrs(length(par.periods))=0;
for i=1:length(pdb)-1
    id1= par.periods>pdb(i)&par.periods<=pdb(i+1)+1;
    smrs(id1)=smris(i);
end
par.smrs=smrs;
userdata.par=par;
set(gcf,'userdata',userdata);

% --- Executes during object creation, after setting all properties.
function smr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in is_figure.
function is_figure_Callback(hObject, eventdata, handles)
% hObject    handle to is_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of is_figure
userdata=get(gcf,'userdata');
par=userdata.par;
par.is_figure=get(handles.is_figure,'value');
userdata.par=par;
set(gcf,'userdata',userdata);

% --- Executes on button press in is_parfor.
function is_parfor_Callback(hObject, eventdata, handles)
% hObject    handle to is_parfor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of is_parfor
userdata=get(gcf,'userdata');
par=userdata.par;
par.is_parfor=get(handles.is_parfor,'value');
userdata.par=par;
set(gcf,'userdata',userdata);

% --- Executes on button press in CIR.
function CIR_Callback(hObject, eventdata, handles)
% hObject    handle to CIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CIR
userdata=get(gcf,'userdata');
par=userdata.par;
par.is_figure=get(handles.is_figure,'value');
userdata.par=par;
set(gcf,'userdata',userdata);



function azmbin_Callback(hObject, eventdata, handles)
% hObject    handle to azmbin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of azmbin as text
%        str2double(get(hObject,'String')) returns contents of azmbin as a double
userdata=get(gcf,'userdata');
par=userdata.par;
par.azmbin=str2double(get(handles.azmbin,'string'));
userdata.par=par;
set(gcf,'userdata',userdata);

% --- Executes during object creation, after setting all properties.
function azmbin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to azmbin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vlim_Callback(hObject, eventdata, handles)
% hObject    handle to vlim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vlim as text
%        str2double(get(hObject,'String')) returns contents of vlim as a double
userdata=get(gcf,'userdata');
par=userdata.par;
par.vlim=str2double(get(handles.vlim,'string'));
userdata.par=par;
set(gcf,'userdata',userdata);

% --- Executes during object creation, after setting all properties.
function vlim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vlim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MaxAzm_Callback(hObject, eventdata, handles)
% hObject    handle to MaxAzm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxAzm as text
%        str2double(get(hObject,'String')) returns contents of MaxAzm as a double
userdata=get(gcf,'userdata');
par=userdata.par;
par.MaxAzm=str2double(get(handles.MaxAzm,'string'));
userdata.par=par;
set(gcf,'userdata',userdata);


% --- Executes during object creation, after setting all properties.
function MaxAzm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxAzm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in onephase.
function onephase_Callback(hObject, eventdata, handles)
% hObject    handle to onephase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of onephase

userdata=get(gcf,'userdata');
par=userdata.par;
par.is_onephase=get(handles.onephase,'value');
userdata.par=par;
set(gcf,'userdata',userdata);


function periods_Callback(hObject, eventdata, handles)
% hObject    handle to periods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of periods as text
%        str2double(get(hObject,'String')) returns contents of periods as a double
userdata=get(gcf,'userdata');
par=userdata.par;


period=get(handles.periods,'string');
period=strsplit(period,':');
if length(period)==1
    try
        period=get(handles.periods,'string');
        period=strsplit(period,',');
        for ipd=1:length(period)
            periods(ipd)=str2double(cell2mat(period(ipd)));
        end
        if ~isempty(find(isnan(periods), 1))
            period=get(handles.periods,'string');
            period=strsplit(period,' ');
            for ipd=1:length(period)
                periods(ipd)=str2double(cell2mat(period(ipd)));
            end
        end
    catch
    end
elseif length(period)==2
    periods=str2double(cell2mat(period(1))):str2double(cell2mat(period(2)));
elseif length(period)==3
    periods=str2double(cell2mat(period(1))):str2double(cell2mat(period(2))):str2double(cell2mat(period(3)));
end

par.periods=periods;
pstr(1:length(par.periods))={nan};
k=0;
for i=par.periods
    k=k+1;
    pstr(k)={i};
end
par.period=par.periods(1);
set(handles.PeriodList,'string',pstr);
set(handles.PeriodList,'value',1);
userdata.par=par;
set(gcf,'userdata',userdata);

% --- Executes during object creation, after setting all properties.
function periods_CreateFcn(hObject, eventdata, handles)
% hObject    handle to periods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in is_bootstrap.
function is_bootstrap_Callback(hObject, eventdata, handles)
% hObject    handle to is_bootstrap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of is_bootstrap
userdata=get(gcf,'userdata');
par=userdata.par;
par.is_bootstrap=get(handles.is_bootstrap,'value');
userdata.par=par;
set(gcf,'userdata',userdata);


function bootTimes_Callback(hObject, eventdata, handles)
% hObject    handle to bootTimes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bootTimes as text
%        str2double(get(hObject,'String')) returns contents of bootTimes as a double
userdata=get(gcf,'userdata');
par=userdata.par;
par.bootTimes=str2double(get(handles.bootTimes,'string'));
userdata.par=par;
set(gcf,'userdata',userdata);

% --- Executes during object creation, after setting all properties.
function bootTimes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bootTimes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function par=getparameters(handles)
userdata=get(gcf,'userdata');
try
par=userdata.par;
catch
end
par.fs='/';
par.indir=get(handles.indir,'string');
par.outdir=get(handles.outdir,'string');
par.outdir2=get(handles.outdir2,'string');
par.bootTimes=str2double(get(handles.bootTimes,'string'));
par.MaxAzm=str2double(get(handles.MaxAzm,'string'));
par.azmbin=str2double(get(handles.azmbin,'string'));
par.vlim=str2double(get(handles.vlim,'string'));
par.is_figure=get(handles.is_figure,'value');
par.is_parfor=get(handles.is_parfor,'value');
par.is_overwrite=get(handles.is_overwrite,'value');
par.is_CIR=get(handles.CIR,'value');
par.is_onephase=get(handles.onephase,'value');
par.is_bootstrap=get(handles.is_bootstrap,'value');
strc=get(handles.strc,'value');
dym=get(handles.dym,'value');

FF=get(handles.faultf,'TooltipString');
FF=strsplit(FF,':');
par.FF=cell2mat(FF(2));

BBF=get(handles.BBF,'TooltipString');
BBF=strsplit(BBF,':');
par.BBF=cell2mat(BBF(2));

if strc==1
    par.vtype='strc';
end
if dym==1
    par.vtype='dym';
end


xylim=get(handles.xylim,'string');
try
    xylim=strsplit(xylim,',');
    Latlim=[str2double(cell2mat(xylim(1))),str2double(cell2mat(xylim(2)))];
    Lonlim=[str2double(cell2mat(xylim(3))),str2double(cell2mat(xylim(4)))];
catch
    xylim=strsplit(xylim,' ');
    Latlim=[str2double(cell2mat(xylim(1))),str2double(cell2mat(xylim(2)))];
    Lonlim=[str2double(cell2mat(xylim(3))),str2double(cell2mat(xylim(4)))];
end
par.Latlim=Latlim;
par.Lonlim=Lonlim;



period=get(handles.periods,'string');
period=strsplit(period,':');
if length(period)==1
    try
        period=get(handles.periods,'string');
        period=strsplit(period,',');
        for ipd=1:length(period)
            periods(ipd)=str2double(cell2mat(period(ipd)));
        end
        if ~isempty(find(isnan(periods), 1))
            period=get(handles.periods,'string');
            period=strsplit(period,' ');
            for ipd=1:length(period)
                periods(ipd)=str2double(cell2mat(period(ipd)));
            end
        end
    catch
    end
elseif length(period)==2
    periods=str2double(cell2mat(period(1))):str2double(cell2mat(period(2)));
elseif length(period)==3
    periods=str2double(cell2mat(period(1))):str2double(cell2mat(period(2))):str2double(cell2mat(period(3)));
end


par.periods=periods;
if ~isfield(par,'period')
    par.period=periods(1);
end

try
    smrTb=get(handles.smrTb,'string');
    smrTb=strsplit(smrTb,',');
    for ism=1:length(smrTb)
        smrTbs(ism)=str2double(cell2mat(smrTb(ism)));
    end
    if ~isempty(find(isnan(smrTbs), 1))
        smrTb=get(handles.smrTb,'string');
        smrTb=strsplit(smrTb,' ');
        for ism=1:length(smrTb)
            smrTbs(ism)=str2double(cell2mat(smrTb(ism)));
        end
    end
catch
end
par.smrTbs=smrTbs;

try
    smr=get(handles.smr,'string');
    smr=strsplit(smr,',');
    for ism=1:length(smr)
        smris(ism)=str2double(cell2mat(smr(ism)));
    end
    if ~isempty(find(isnan(smris), 1))
        smr=get(handles.smr,'string');
        smr=strsplit(smr,' ');
        for ism=1:length(smr)
            smris(ism)=str2double(cell2mat(smr(ism)));
        end
    end
catch
end
pdb=[min(par.periods)-1,par.smrTbs,max(par.periods)+1];
smrs(length(par.periods))=0;
for i=1:length(pdb)-1
    id1= par.periods>pdb(i)&par.periods<=pdb(i+1)+1;
    smrs(id1)=smris(i);
end
par.smrs=smrs;

#load coast
#par.long=long;
#par.lat=lat;

load coastlines
par.long=coastlon;
par.lat=coastlat;

% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
