function varargout = STdispersion(varargin)
% STDISPERSION MATLAB code for STdispersion.fig
%      STDISPERSION, by itself, creates a new STDISPERSION or raises the existing
%      singleton*.
%
%      H = STDISPERSION returns the handle to a new STDISPERSION or the handle to
%      the existing singleton*.
%
%      STDISPERSION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STDISPERSION.M with the given input arguments.
%
%      STDISPERSION('Property','Value',...) creates a new STDISPERSION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before STdispersion_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to STdispersion_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help STdispersion

% Last Modified by GUIDE v2.5 13-Aug-2023 13:33:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @STdispersion_OpeningFcn, ...
    'gui_OutputFcn',  @STdispersion_OutputFcn, ...
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


% --- Executes just before STdispersion is made visible.
function STdispersion_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to STdispersion (see VARARGIN)

% Choose default command line output for STdispersion
handles.output = hObject;
indir='E:/ETP/AutoWG_gridV5/aniso/station/';
outdir='E:/ETP/AutoWG_gridV5/aniso/station/Pick/';
par.isoverwrite=0;
stfs=dir([indir 'bk*.mat']);



par.indir=indir;
par.outdir=outdir;
par.dvm=str2double(get(handles.paradvm,'string'));
par.dviso=str2double(get(handles.paradviso,'string'));
par.dvtance=str2double(get(handles.paradvtance,'string'));
par.dvazm=str2double(get(handles.paradvazm,'string'));

par.stfs=stfs;
for sti=1:length(stfs)
    stfi=stfs(sti).name;
    strst=strsplit(stfi,'_');
    xi=str2double(cell2mat(strst(2)));
    yi=str2double(cell2mat(strsplit(cell2mat(strst(3)),'.mat')));    
    par.stloc(sti,:)=[xi,yi];
    par.stname(sti)={stfi};    
    try
        load([outdir '/' stfi])        
    catch
        load([indir '/' stfi])
    end    
    for pi=1:length(wgst.periods)
        vi=wgst.v(:,pi);
        par.viso(sti,pi)=median(vi(~isnan(vi)));
    end       
end


for i=1:length(par.viso(1,:))
    visoi=par.viso(:,i);
    par.visom(i)=median(visoi(~isnan(visoi)));
end

set(handles.stlist,'string',par.stname);
set(handles.stlist,'value',1);
par.handles=handles;
par.g=1;
par.h=1;
par.outfile=[par.outdir '/' par.stfs(1).name];
try
load(par.outfile)
catch
load([par.indir par.stfs(1).name])
end

load coast
par.long=long;
par.lat=lat;

stj=wgst;
azmLim=str2double(get(handles.azmLim,'string'));
nbgrd=str2double(get(handles.nbgrd,'string'));
hg=findnearDisp(par,nbgrd,azmLim,stj);
par.hg=hg;

strazm=round(azimuth(stj.evla,stj.evlo,stj.st(1),stj.st(2)));
set(handles.evlist,'string',strazm);
set(handles.evlist,'value',1);
% IMGstaDisp(par);
azmLim=str2double(get(handles.azmLim,'string'));
nbgrd=str2double(get(handles.nbgrd,'string'));
hg=findnearDisp(par,nbgrd,azmLim,stj);
par.hg=hg;
% PlotStationLocation(handles,par)
% plotazm(par,stj)

userdata.par=par;
set(gcf,'userdata',userdata);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes STdispersion wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = STdispersion_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in stlist.
function stlist_Callback(hObject, eventdata, handles)
% hObject    handle to stlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns stlist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from stlist
userdata=get(gcf,'userdata');
par = userdata.par;
par.g=1;
par.h=get(handles.stlist,'value');
par.dvm=str2double(get(handles.paradvm,'string'));
par.dviso=str2double(get(handles.paradviso,'string'));
par.dvtance=str2double(get(handles.paradvtance,'string'));
par.dvazm=str2double(get(handles.paradvazm,'string'));
par.outfile=[par.outdir '/' par.stfs(par.h).name];
try
load(par.outfile)
catch
load([par.indir '/' par.stfs(par.h).name])
end

stj=wgst;
azmLim=str2double(get(handles.azmLim,'string'));
nbgrd=str2double(get(handles.nbgrd,'string'));
hg=findnearDisp(par,nbgrd,azmLim,stj);
par.hg=hg;

strazm=round(azimuth(stj.evla,stj.evlo,stj.st(1),stj.st(2)));
set(handles.evlist,'string',strazm);
set(handles.evlist,'value',1);
IMGstaDisp(par);
subplot(handles.stmap);
hold on
PlotStationLocation(handles,par)
plotazm(par,stj)
% plot(par.stloc(:,2), par.stloc(:,1),'ro');
hold off
userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function stlist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in evlist.
function evlist_Callback(hObject, eventdata, handles)
% hObject    handle to evlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns evlist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from evlist

userdata=get(gcf,'userdata');
par = userdata.par;
par.dvm=str2double(get(handles.paradvm,'string'));
par.dviso=str2double(get(handles.paradviso,'string'));
par.dvtance=str2double(get(handles.paradvtance,'string'));
par.dvazm=str2double(get(handles.paradvazm,'string'));
par.g=get(handles.evlist,'value');
par=plotdisp(par);
plotazm(par,par.hg.stj)
userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function evlist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to evlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RmEv.
function RmEv_Callback(hObject, eventdata, handles)
% hObject    handle to RmEv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par = userdata.par;
par.dvm=str2double(get(handles.paradvm,'string'));
par.dviso=str2double(get(handles.paradviso,'string'));
par.dvtance=str2double(get(handles.paradvtance,'string'));
par.dvazm=str2double(get(handles.paradvazm,'string'));
par.g=get(handles.evlist,'value');
wgst=par.hg.stj;
wgst.evla(par.g)=[];
wgst.evlo(par.g)=[];
wgst.v(par.g,:)=[];
wgst.va(par.g,:)=[];
wgst.vcr(par.g,:)=[];
wgst.vcra(par.g,:)=[];
wgst.dvcr(par.g,:)=[];
wgst.weight(par.g,:)=[];
wgst.azm(par.g,:)=[];
wgst.Bx(par.g,:)=[];
wgst.By(par.g,:)=[];
save(par.outfile,'wgst')
par.hg.stj=wgst;
if par.g>length(wgst.evla)
    par.g=length(wgst.evla);
end

azm=round(azimuth(wgst.evla,wgst.evlo,wgst.st(1),wgst.st(2)));
par.azm=azm;
set(handles.evlist,'string',azm);


set(handles.evlist,'value',par.g);
par=plotdisp(par);
PlotStationLocation(handles,par)
plotazm(par,wgst)

userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);

% --- Executes on button press in RmEv.
function RmSt_Callback(hObject, eventdata, handles)
% hObject    handle to RmEv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par = userdata.par;
par.dvm=str2double(get(handles.paradvm,'string'));
par.dviso=str2double(get(handles.paradviso,'string'));
par.dvtance=str2double(get(handles.paradvtance,'string'));
par.dvazm=str2double(get(handles.paradvazm,'string'));
par.h=get(handles.stlist,'value');

if ~exist([par.indir '/bad/'],'dir')
    mkdir([par.indir '/bad/'])
end

movefile([par.indir par.stfs(par.h).name],[par.indir '/bad/'])
par.stfs(par.h)=[];
par.stname(par.h)=[];
if par.h>length(par.stfs)
   par.h=length(par.stfs);
else
    par.h=par.h;
end

set(handles.stlist,'string',par.stname);
set(handles.stlist,'value',par.h);
par.handles=handles;

par.outfile=[par.outdir '/' par.stfs(par.h).name];
try
load(par.outfile)
catch
load([par.indir par.stfs(1).name])
end

load coast
par.long=long;
par.lat=lat;

stj=wgst;
azmLim=str2double(get(handles.azmLim,'string'));
nbgrd=str2double(get(handles.nbgrd,'string'));
hg=findnearDisp(par,nbgrd,azmLim,stj);
par.hg=hg;

strazm=round(azimuth(stj.evla,stj.evlo,stj.st(1),stj.st(2)));
set(handles.evlist,'string',strazm);
set(handles.evlist,'value',1);
IMGstaDisp(par);
azmLim=str2double(get(handles.azmLim,'string'));
nbgrd=str2double(get(handles.nbgrd,'string'));
hg=findnearDisp(par,nbgrd,azmLim,stj);
par.hg=hg;
PlotStationLocation(handles,par)
plotazm(par)

userdata.par=par;
set(gcf,'userdata',userdata);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in SAVE.
function SAVE_Callback(hObject, eventdata, handles)
% hObject    handle to SAVE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par=userdata.par;
disp('done')


% --- Executes on button press in Dr.
function Dr_Callback(hObject, eventdata, handles)
% hObject    handle to Dr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par=userdata.par;
par.dvm=str2double(get(handles.paradvm,'string'));
par.dviso=str2double(get(handles.paradviso,'string'));
par.dvtance=str2double(get(handles.paradvtance,'string'));
par.dvazm=str2double(get(handles.paradvazm,'string'));
wgst=par.hg.stj;
par.PdlY=[find(par.deleteID) par.PdlY];
par.PdlY=unique(par.PdlY);
wgst.v(par.g,par.PdlY)=nan;
wgst.va(par.g,par.PdlY)=nan;
wgst.vcr(par.g,par.PdlY)=nan;
wgst.vcra(par.g,par.PdlY)=nan;
wgst.dvcr(par.g,par.PdlY)=nan;
wgst.weight(par.g,par.PdlY)=nan;
wgst.azm(par.g,par.PdlY)=nan;
wgst.Bx(par.g,par.PdlY)=nan;
wgst.By(par.g,par.PdlY)=nan;
par.hg.stj=wgst;
save(par.outfile,'wgst')
par.PdlX=[];
par.PdlY=[];
par=plotdisp(par);
% Update handles structure
userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);

% --- Executes on button press in Dg.
function Dg_Callback(hObject, eventdata, handles)
% hObject    handle to Dg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par=userdata.par;
par.dvm=str2double(get(handles.paradvm,'string'));
par.dviso=str2double(get(handles.paradviso,'string'));
par.dvtance=str2double(get(handles.paradvtance,'string'));
par.dvazm=str2double(get(handles.paradvazm,'string'));
par.PdlX=[find(par.deleteID) par.PdlX];
par.PdlX=unique(par.PdlX);
wgst=par.hg.stj;
wgst.v(par.g,par.PdlX)=nan;
wgst.va(par.g,par.PdlX)=nan;
wgst.vcr(par.g,par.PdlX)=nan;
wgst.vcra(par.g,par.PdlX)=nan;
wgst.dvcr(par.g,par.PdlX)=nan;
wgst.weight(par.g,par.PdlX)=nan;
wgst.azm(par.g,par.PdlY)=nan;
wgst.Bx(par.g,par.PdlY)=nan;
wgst.By(par.g,par.PdlY)=nan;
par.hg.stj=wgst;
save(par.outfile,'wgst')
par.PdlX=[];
par.PdlY=[];
par=plotdisp(par);
% Update handles structure
userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);

% --- Executes on button press in Lev.
function Lev_Callback(hObject, eventdata, handles)
% hObject    handle to Lev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par = userdata.par;
par.dvm=str2double(get(handles.paradvm,'string'));
par.dviso=str2double(get(handles.paradviso,'string'));
par.dvtance=str2double(get(handles.paradvtance,'string'));
par.dvazm=str2double(get(handles.paradvazm,'string'));
par.g=get(handles.evlist,'value');

if par.g>1
    par.g=par.g-1;
end
set(handles.evlist,'value',par.g);

par=plotdisp(par);
plotazm(par,par.hg.stj)

userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);


% --- Executes on button press in NextEv.
function NextEv_Callback(hObject, eventdata, handles)
% hObject    handle to NextEv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par = userdata.par;
par.dvm=str2double(get(handles.paradvm,'string'));
par.dviso=str2double(get(handles.paradviso,'string'));
par.dvtance=str2double(get(handles.paradvtance,'string'));
par.dvazm=str2double(get(handles.paradvazm,'string'));
par.g=get(handles.evlist,'value');
stj=par.hg.stj;
if par.g<length(stj.evlo)
    par.g=par.g+1;
end

par=plotdisp(par);
plotazm(par,par.hg.stj)
set(handles.evlist,'value',par.g);
userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);


% --- Executes on button press in Pick.
function Pick_Callback(hObject, eventdata, handles)
% hObject    handle to Pick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par=userdata.par;
par.dvm=str2double(get(handles.paradvm,'string'));
par.dviso=str2double(get(handles.paradviso,'string'));
par.dvtance=str2double(get(handles.paradvtance,'string'));
par.dvazm=str2double(get(handles.paradvazm,'string'));
par=plotdisp(par);
hold(handles.figdsp1,'on')
par.dlX=[];
par.dlY=[];
v=par.hg.stj.v(par.g,:);
periods=par.hg.stj.periods;
[y,x] = getpts(handles.figdsp1);
% [ys,xs] = getpts(handles.figdsp1);
% xidx=1:length(xs);
% id1=mod(xidx,2)==1;
% id2=mod(xidx,2)==0;

% % [y,x] = ginput(2);
% for i=1:length(find(mod(xidx,2)==1))
deleteXi=find(v>min(x)&v<max(x)&periods>min(y)&periods<max(y));
deleteYi=find(v>min(x)&v<max(x)&periods>min(y)&periods<max(y));
% end
scatter(periods(deleteXi),v(deleteXi),25,'fill','g')
hold(handles.figdsp1,'on')
scatter(periods(deleteYi),v(deleteYi),20,'fill','r')
hold(handles.figdsp1,'off')
par.PdlX=deleteXi;
par.PdlY=deleteYi;
userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);

% --- Executes on button press in Auto.
function Auto_Callback(hObject, eventdata, handles)
% hObject    handle to Auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userdata=get(gcf,'userdata');
par = userdata.par;
par.handles=[];
par.dvm=str2double(get(handles.paradvm,'string'));
par.dviso=str2double(get(handles.paradviso,'string'));
par.dvtance=str2double(get(handles.paradvtance,'string'));
par.dvazm=str2double(get(handles.paradvazm,'string'));
azmLim=str2double(get(handles.azmLim,'string'));
nbgrd=str2double(get(handles.nbgrd,'string'));
loop=get(handles.Loops,'string');
loop=str2double(loop);
paramat='paraDispPick.mat';
save(paramat,'par','azmLim','nbgrd','loop')

disp('Auto Pickign...')
for loops=1:loop
    parfor h=1:length(par.stfs)
        autoPick(paramat,h)
    end  
end
disp('Picking done')
delete(paramat)



function autoPick(paramat,h)
load(paramat)
par.h=h;
if mod(h/length(par.stfs)*100,10)==0
    disp([num2str(h/length(par.stfs)*100),'%'])
end


infile=[par.indir '/' par.stfs(h).name];
outfile=[par.outdir '/' par.stfs(h).name];
if exist(outfile,'file')&&par.isoverwrite==0;
    disp(['overwrite=0, skip:  ' par.stfs(h).name])
    return
end

try
    load(outfile)
catch
    load(infile)
end
stj=wgst;
igk=0;
for g=1:length(stj.evla)
    par.g=g;
    vr=stj.v(g,:);
    periods=stj.periods;
    v=vr;
    prd=stj.periods;
    
    if isempty(v(~isnan(v)))
        igk=igk+1;
        badEVi(igk)=g;
        continue
    end
    
    idv1=abs(v-par.visom);
    bad1= idv1>par.dvm|isnan(v);
    %%
    hg=findnearDisp(par,nbgrd,azmLim,stj);
    par.hg=hg;
    hi=hg.hi;
    gii=hg.gii;
    for i=1:length(hi)
        idg=cell2mat(gii(i));
        sti=hg.sti(i);
        vi(1:length(idg),:)=sti.v(idg,:);
        vis=sti.v;
        for ip=1:length(vi(1,:))
            vmip=vi(:,ip);
            vmips=vis(:,ip);
            vmip(isnan(vmip))=[];
            vmips(isnan(vmips))=[];
            if ~isempty(vmip)
                vmi(i,ip)=median(vmip);
            else
                vmi(i,ip)=nan;
            end
            
            if ~isempty(vmips)
                vmis(i,ip)=median(vmips);
            else
                vmis(i,ip)=nan;
            end
            
        end
    end
    vmis=[vmis;par.viso(h,:)];
    
    for ip=1:length(stj.periods)
        vmip=vmi(:,ip);
        vmip(isnan(vmip))=[];
        vm(ip)=median(vmip);
        vmips=vmis(:,ip);
        vmips(isnan(vmips))=[];
        vmiso(ip)=median(vmips);
    end
    
    prdvm=stj.periods;
    prdvm(isnan(vm))=[];
    vm(isnan(vm))=[];
    
    if length(vm)>=10
        pvm = polyfit(prdvm,vm,3);
        vvm=polyval(pvm,prd);
    else
        vvm=vm;
    end
    
    
    idv4=abs(v-vvm);
    bad4= idv4>par.dvazm|isnan(v);
    
    
    idv2=abs(v-vmiso);
    bad2= idv2>par.dviso|isnan(v);
    
    
    gid=~bad1&~bad2&~bad4;
    
    if length(v(gid))>=10
        pv = polyfit(prd(gid),v(gid),3);
        vv=polyval(pv,prd);
    elseif isempty(v(gid))
        continue
    else
        vv=v;
    end
    
    
    idv3=abs(v-vv);
    bad3= idv3>par.dvtance|isnan(v);
    
    badid = bad1 | bad2 | bad3 | bad4;
    
    stj.v(g,badid)=nan;
    stj.va(g,badid)=nan;
    stj.vcr(g,badid)=nan;
    stj.vcra(g,badid)=nan;
    stj.dvcr(g,badid)=nan;
    stj.weight(g,badid)=nan;
    stj.azm(g,badid)=nan;
    stj.Bx(g,badid)=nan;
    stj.By(g,badid)=nan;
    
    if isempty(stj.v(~isnan(stj.v)))
        igk=igk+1;
        badEVi(igk)=g;
    end
    
end

if igk>0
    stj.evla(badEVi)=[];
    stj.evlo(badEVi)=[];
    stj.v(badEVi,:)=[];
    stj.va(badEVi,:)=[];
    stj.vcr(badEVi,:)=[];
    stj.vcra(badEVi,:)=[];
    stj.dvcr(badEVi,:)=[];
    stj.weight(badEVi,:)=[];
    stj.azm(badEVi,:)=[];
    stj.Bx(badEVi,:)=[];
    stj.By(badEVi,:)=[];
end
wgst=stj;
save(outfile,'wgst')






% --- Executes on button press in stnext.
function stnext_Callback(hObject, eventdata, handles)
% hObject    handle to stnext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par = userdata.par;
par.dvm=str2double(get(handles.paradvm,'string'));
par.dviso=str2double(get(handles.paradviso,'string'));
par.dvtance=str2double(get(handles.paradvtance,'string'));
par.dvazm=str2double(get(handles.paradvazm,'string'));
par.h=get(handles.stlist,'value');
par.g=1;
if par.h<length(par.stfs)
    par.h=par.h+1;
end
par.outfile=[par.outdir '/' par.stfs(par.h).name];
try
load(par.outfile)
catch
load([par.indir par.stfs(par.h).name])
end

stj=wgst;
azmLim=str2double(get(handles.azmLim,'string'));
nbgrd=str2double(get(handles.nbgrd,'string'));
hg=findnearDisp(par,nbgrd,azmLim,stj);
par.hg=hg;

strazm=round(azimuth(stj.evla,stj.evlo,stj.st(1),stj.st(2)));
set(handles.evlist,'string',strazm);
set(handles.evlist,'value',1);
set(handles.stlist,'value',par.h);

IMGstaDisp(par);
subplot(handles.stmap);
hold on
PlotStationLocation(handles,par)
plotazm(par,stj)
% plot(par.stloc(:,2), par.stloc(:,1),'ro');
hold off
userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);


% --- Executes on button press in stLast.
function stLast_Callback(hObject, eventdata, handles)
% hObject    handle to stLast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par = userdata.par;
par.dvm=str2double(get(handles.paradvm,'string'));
par.dviso=str2double(get(handles.paradviso,'string'));
par.dvtance=str2double(get(handles.paradvtance,'string'));
par.dvazm=str2double(get(handles.paradvazm,'string'));
par.h=get(handles.stlist,'value');
par.g=1;
if par.h>1
    par.h=par.h-1;
end
par.outfile=[par.outdir '/' par.stfs(par.h).name];
try
load(par.outfile)
catch
load([par.indir par.stfs(par.h).name])
end

stj=wgst;
azmLim=str2double(get(handles.azmLim,'string'));
nbgrd=str2double(get(handles.nbgrd,'string'));
hg=findnearDisp(par,nbgrd,azmLim,stj);
par.hg=hg;

strazm=round(azimuth(stj.evla,stj.evlo,stj.st(1),stj.st(2)));
set(handles.evlist,'string',strazm);
set(handles.evlist,'value',1);
set(handles.stlist,'value',par.h);

IMGstaDisp(par);
subplot(handles.stmap);
hold on
PlotStationLocation(handles,par)
plotazm(par,stj)
% plot(par.stloc(:,2), par.stloc(:,1),'ro');
hold off
userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);

function PlotStationLocation(handles,par)
h1 = handles.stmap;
hold(h1, 'off');
g_allsta.lat=par.stloc(:,2);
g_allsta.lon=par.stloc(:,1);
plot(h1, g_allsta.lon, g_allsta.lat, '.','color',[0.6 0.6 0.6]);
hold(h1, 'on');
hg=par.hg;
hi=hg.hi;
plot(h1, g_allsta.lon(hi), g_allsta.lat(hi), 'r.');
scatter(h1, g_allsta.lon(hg.h), g_allsta.lat(hg.h), 'ob')
axis equal
% set(h1,'Visible','off')



function IMGstaDisp(par)
handles=par.handles;
stj=par.hg.stj;
f1=handles.figdsp1;
hold(f1, 'off');

for i=1:length(stj.evla)
    periods=stj.periods;
    v=stj.v(i,:);
    periods(isnan(v))=[];
    v(isnan(v))=[];
    
    if ~isempty(v)
        vv=v;
        azm=stj.azm(i,1);        
        
        plot(f1,periods,vv,'b.')
        
        xlim(f1,[min(stj.periods)-5,max(stj.periods)+5])
        hold(f1, 'on');
    else
        plot(f1,[min(stj.periods)-5,max(stj.periods)+5],[4.5,4.5],'k')
    end
end
ylim(f1,[2,5])
hold(f1, 'off');





function par=plotdisp(par)
handles=par.handles;
st=par.hg.stj;
periods=st.periods;
v=st.v(par.g,:);
dv3=abs(v-par.visom);

%%
azmLim=str2double(get(handles.azmLim,'string'));
nbgrd=str2double(get(handles.nbgrd,'string'));
hg=findnearDisp(par,nbgrd,azmLim,st);
par.hg=hg;
hi=hg.hi;
gii=hg.gii;
f2=handles.figdsp1;
for i=1:length(hi)    
    idg=cell2mat(gii(i));
    sti=hg.sti(i);
    vi(1:length(idg),:)=sti.v(idg,:);
    vis=sti.v;
    for ip=1:length(vi(1,:))
        vmip=vi(:,ip);
        vmip(isnan(vmip))=[];
        if ~isempty(vmip)
            vmi(i,ip)=median(vmip);
        else
            vmi(i,ip)=nan;
        end
        
        vmips=vis(:,ip);
        vmips(isnan(vmips))=[];
        if ~isempty(vmips)
            vmis(i,ip)=median(vmips);
        else
            vmis(i,ip)=nan;
        end
        
    end    
    plot(f2,st.periods,vmi(i,:),'g')
    hold(f2,'on')

end
%%
vmis=[vmis;par.viso(par.h,:)];
for ip=1:length(st.periods)
    vmipi=vmi(:,ip);    
    vmipi(isnan(vmipi))=[];
    vm(ip)=median(vmipi);
    vmipis=vmis(:,ip);    
    vmipis(isnan(vmipis))=[];
    vms(ip)=median(vmipis);
end
prdvm=st.periods;
prdvm(isnan(vm))=[];
vm(isnan(vm))=[];
vms(isnan(vms))=[];

dv2=abs(v-vms);

pvm = polyfit(prdvm,vm,3);
vvm=polyval(pvm,periods);
plot(f2,periods,vm,'.k','LineWidth',16)
dv4=abs(v-vm);   

f2=handles.figdsp1;
plot(f2,st.periods,par.viso(par.h,:),'r','LineWidth',3)
plot(f2,st.periods,par.visom,'y','LineWidth',3)
plot(f2,st.periods,par.viso(par.h,:)-par.dvm,'b')
plot(f2,st.periods,par.viso(par.h,:)+par.dvm,'b')
scatter(f2,periods,v,40,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0 .7 .8])
idv=dv2<=par.dviso&dv3<=par.dvm&dv4<=par.dvazm;
if length(v)>5
    pv = polyfit(periods(idv),v(idv),3);
    vv=polyval(pv,periods,3);
else
    vv = v;
end
plot(f2,periods(~isnan(v)),vv(~isnan(v)),'k','LineWidth',1);
dv1=abs(v-vv);


scatter(f2,periods(dv1>par.dvtance),v(dv1>par.dvtance),20,'MarkerEdgeColor','w',...
    'MarkerFaceColor','R')
scatter(f2,periods(dv2>par.dviso),v(dv2>par.dviso),20,'MarkerEdgeColor','w',...
    'MarkerFaceColor','R')
scatter(f2,periods(dv3>par.dvm),v(dv3>par.dvm),20,'MarkerEdgeColor','w',...
    'MarkerFaceColor','R')
scatter(f2,periods(dv4>par.dvazm),v(dv4>par.dvazm),20,'MarkerEdgeColor','w',...
    'MarkerFaceColor','R')

par.deleteID=dv1>par.dvtance|dv2>par.dviso|dv3>par.dvm|dv4>par.dvazm;

xlim(f2,[min(st.periods)-5,max(st.periods)+5])
ylim(f2,[2,5])
hold(f2,'off')
par.Pvd=periods(dv1>par.dvtance);
par.Pvd2=periods(dv2>par.dvm);


function hg=findnearDisp(par,KmLim1,azmLim,stj)
h=par.h;
g=par.g;
hg.stj=stj;
azm=azimuth(stj.evla,stj.evlo,stj.st(1),stj.st(2));
azm0=azm(g);
par.azm=azm;
dazm=azm-azm0;
while ~isempty(find(abs(dazm)>90, 1))
dazm(dazm>90)=dazm(dazm>90)-180;   
dazm(dazm<-90)=dazm(dazm<-90)+180;   
end
hg.gi=find(dazm>-azmLim&dazm<azmLim);
hg.h=h;


dist=distance(stj.st(1),stj.st(2),par.stloc(:,2),par.stloc(:,1));
hg.hi=find(dist>-KmLim1/deg2km(1)&dist<KmLim1/deg2km(1));

for i=1:length(hg.hi)
    index=hg.hi(i);
    stf=par.stfs(index).name;
    
    try
        load([par.indir '/Pick_' stf])
    catch
        load([par.indir '/' stf])
    end
    
    sti=wgst;
    hg.sti(i)=wgst;
    azm=azimuth(sti.evla,sti.evlo,sti.st(1),sti.st(2));
    par.azm=azm;
    dazm=azm-azm0;
    while ~isempty(find(abs(dazm)>90, 1))
        dazm(dazm>90)=dazm(dazm>90)-180;
        dazm(dazm<-90)=dazm(dazm<-90)+180;
    end
    hg.gii(i,1)={find(dazm>-azmLim&dazm<azmLim)};
end



function plotazm(par,stj)
handles=par.handles;
par.g=get(handles.evlist,'value');
par.h=get(handles.stlist,'value');
azmLim=str2double(get(handles.azmLim,'string'));
nbgrd=str2double(get(handles.nbgrd,'string'));
hg=findnearDisp(par,nbgrd,azmLim,stj);
par.hg=hg;
evlo=stj.evlo(par.g);
evla=stj.evla(par.g);
stla=stj.st(1);
stlo=stj.st(2);
af=handles.azmap;
%%world map

plot(af,par.long,par.lat,'Color',[150 150 150]/255)
hold(af,'on')
plot(af,stj.evlo,stj.evla,'.r')
scatter(af,stlo,stla,40,'^','MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1)
xlim(af,[-180 180])
ylim(af,[-90 90])
hi=hg.hi;
gii=hg.gii;

for i=1:length(hi)
    stgii=hg.sti(i);    
    plot(af,stgii.evlo(cell2mat(gii(i))),stgii.evla(cell2mat(gii(i))),'.g')    
end

stgi=hg.stj; 
gi=hg.gi;
plot(af,stgi.evlo(gi),stgi.evla(gi),'.g')
scatter(af,evlo,evla,'b*')
hold(af,'off')
% set(af,'Visible','off')



function Loops_Callback(hObject, eventdata, handles)
% hObject    handle to Loops (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Loops as text
%        str2double(get(hObject,'String')) returns contents of Loops as a double


% --- Executes during object creation, after setting all properties.
function Loops_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Loops (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function azmLim_Callback(hObject, eventdata, handles)
% hObject    handle to azmLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of azmLim as text
%        str2double(get(hObject,'String')) returns contents of azmLim as a double


% --- Executes during object creation, after setting all properties.
function azmLim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to azmLim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nbgrd_Callback(hObject, eventdata, handles)
% hObject    handle to nbgrd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nbgrd as text
%        str2double(get(hObject,'String')) returns contents of nbgrd as a double


% --- Executes during object creation, after setting all properties.
function nbgrd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nbgrd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function paradvm_Callback(hObject, eventdata, handles)
% hObject    handle to paradvm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of paradvm as text
%        str2double(get(hObject,'String')) returns contents of paradvm as a double


% --- Executes during object creation, after setting all properties.
function paradvm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to paradvm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function paradviso_Callback(hObject, eventdata, handles)
% hObject    handle to paradviso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of paradviso as text
%        str2double(get(hObject,'String')) returns contents of paradviso as a double


% --- Executes during object creation, after setting all properties.
function paradviso_CreateFcn(hObject, eventdata, handles)
% hObject    handle to paradviso (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function paradvtance_Callback(hObject, eventdata, handles)
% hObject    handle to paradvtance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of paradvtance as text
%        str2double(get(hObject,'String')) returns contents of paradvtance as a double


% --- Executes during object creation, after setting all properties.
function paradvtance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to paradvtance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function paradvazm_Callback(hObject, eventdata, handles)
% hObject    handle to paradvazm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of paradvazm as text
%        str2double(get(hObject,'String')) returns contents of paradvazm as a double


% --- Executes during object creation, after setting all properties.
function paradvazm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to paradvazm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if isempty(eventdata.Modifier)%判断是否有快捷键
    ctrl='';%无则赋值为空
else
    ctrl=eventdata.Modifier{1};%有，则取出辅助键
end
key=eventdata.Key;
switch key
    case 'alt'
        Pick_Callback(handles.Pick,[],handles)
    case 'downarrow'
        stnext_Callback(handles.stnext,[],handles)
    case 'uparrow'
        stLast_Callback(handles.stLast,[],handles)
    case 'r'
        Dr_Callback(handles.Dr,[],handles)
    case 'g'
        Dg_Callback(handles.Dg,[],handles)
    case 'space'
        NextEv_Callback(hObject, eventdata, handles)        
    otherwise
end


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par =userdata.par;
switch (get(gcbf,'SelectionType'))
    case 'alt'
        Pick_Callback(handles.Pick,[],handles)
    case 'normal'
%         par=plotdisp(par);
    otherwise
end

% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)
count = eventdata.VerticalScrollCount;
if count==-1
    Lev_Callback(hObject, eventdata, handles)
elseif count==1
    NextEv_Callback(hObject, eventdata, handles)
end


% --- Executes on mouse press over axes background.
function figdsp1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figdsp1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
test=1;
