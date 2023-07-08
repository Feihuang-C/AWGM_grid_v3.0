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

% Last Modified by GUIDE v2.5 07-Sep-2021 20:32:15

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
% temp=load('D:\dyx\cfh\AutoGradANA_Grid3\aniso\wga.R.BHZ-dvazm0.3.mat');22518testwga.R.BHZ
% par.outfile='D:\dyx\cfh\AutoGradANA_Grid3\aniso\wga.R.BHZ-dvazm0.15.mat';
temp=load('./aniso/22518testwga.R.BHZ.mat');
par.outfile='./aniso/testplot.mat';
% temp=load('D:\dyx\cfh\AutoGradANA_Grid3\aniso\0706teststdNETibetAniso2180ban10vlim0.35.R.BHZ.mat');
% par.outfile='D:\dyx\cfh\AutoGradANA_Grid3\aniso\STd22518testwga.R.BHZ.mat';
par.periods=10:80;

% temp=load('J:\wga.R.BHZ.mat');
% par.outfile='J:\Qwga.R.BHZ.mat';
% par.periods=10:120;

par.dvm=str2double(get(handles.paradvm,'string'));
par.dviso=str2double(get(handles.paradviso,'string'));
par.dvtance=str2double(get(handles.paradvtance,'string'));
par.dvazm=str2double(get(handles.paradvazm,'string'));
st=temp.stk;
k=0;
for sti=1:length(st)
    if isempty(st(sti).evla)
        k=k+1;
        emti(k)=sti;
    end
end
if k>0
    st(emti)=[];
end

for sti=1:length(st)
    stname(sti)=st(sti).stn;
    par.stloc(sti,:)=[st(sti).stla,st(sti).stlo];
    for i=1:length(par.periods);
        visom=st(sti).v0(:,i);
        visom(isnan(visom))=[];
        vm(i)=median(visom);
    end
    prd=par.periods(~isnan(vm));
    vm(isnan(vm))=[];
    pv = polyfit(prd,vm,3);
    par.viso(sti,:)=polyval(pv,par.periods);
end

par.stname=stname;
set(handles.stlist,'string',stname);
set(handles.stlist,'value',1);

for i=1:length(par.periods);
    visom=par.viso(:,i);
    visom(isnan(visom))=[];
    vmall(i)=median(visom);
end
prd=par.periods(~isnan(vmall));
vmall(isnan(vmall))=[];
pvmall = polyfit(prd,vmall,3);
par.vmall=polyval(pvmall,par.periods);


par.handles=handles;
par.st=st;
par.g=1;
par.h=1;

stj=par.st(par.h);
for evi=1:length(stj.evla)
    azm(evi)=azimuth0(median(par.stloc(:,1)),median(par.stloc(:,2)),stj.evla(evi),stj.evlo(evi));
end
par.azm=azm;
for evi=1:length(azm)
    strazm(evi)={['azm' num2str2(azm(evi),5,1)]};
end

set(handles.evlist,'string',strazm);
set(handles.evlist,'value',1);
par=IMGstaDisp(par);
stmp=handles.stmap;
hold(stmp,'on')
stlat=par.st(par.h).stla;
stlon=par.st(par.h).stlo;
plot(stmp,stlon,stlat,'ro');
hold(stmp,'off')

azmLim=str2double(get(handles.azmLim,'string'));
nbgrd=str2double(get(handles.nbgrd,'string'));
hg=findnearDisp(par,nbgrd,azmLim);
par.hg=hg;
PlotStationLocation(handles,par)
plotazm(par)

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
par.h=get(handles.stlist,'value');
par.dvm=str2double(get(handles.paradvm,'string'));
par.dviso=str2double(get(handles.paradviso,'string'));
par.dvtance=str2double(get(handles.paradvtance,'string'));;
par.dvazm=str2double(get(handles.paradvazm,'string'));

stj=par.st(par.h);
for evi=1:length(stj.evla)
    azm(evi)=azimuth0(stj.stla,stj.stlo,stj.evla(evi),stj.evlo(evi));
end
par.azm=azm;
for evi=1:length(azm)
    strazm(evi)={['bazm' num2str2(azm(evi),5,1)]};
end
set(handles.evlist,'string',strazm);
set(handles.evlist,'value',1);
par=IMGstaDisp(par);
subplot(handles.stmap);
hold on

azmLim=str2double(get(handles.azmLim,'string'));
nbgrd=str2double(get(handles.nbgrd,'string'));
hg=findnearDisp(par,nbgrd,azmLim);
par.hg=hg;

PlotStationLocation(handles,par)
plotazm(par)
stlat=par.st(par.h).stla;
stlon=par.st(par.h).stlo;
plot(stlon,stlat,'ro');
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
par.dvtance=str2double(get(handles.paradvtance,'string'));;
par.dvazm=str2double(get(handles.paradvazm,'string'));
par.g=get(handles.evlist,'value');
par=plotdisp(par);
plotazm(par)
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
st=par.st;
st(par.h).evla(par.g)=[];
st(par.h).evlo(par.g)=[];
st(par.h).dpth(par.g)=[];
st(par.h).v0(par.g,:)=[];
st(par.h).dv(par.g,:)=[];
st(par.h).a0(par.g,:)=[];
st(par.h).da(par.g,:)=[];
st(par.h).r0(par.g,:)=[];
st(par.h).dr(par.g,:)=[];
st(par.h).g0(par.g,:)=[];
st(par.h).dg(par.g,:)=[];
st(par.h).vg(par.g,:)=[];
st(par.h).azmo(par.g,:)=[];
st(par.h).periods(par.g,:)=[];
st(par.h).nsti(par.g,:)=[];
par.st=st;
stj=st(par.h);
if par.g>length(stj.v0(:,1))
    par.g=length(stj.v0(:,1));
end

for evi=1:length(stj.evla)
    azm(evi)=azimuth0(st(par.h).stla,st(par.h).stlo,st(par.h).evla(evi),st(par.h).evlo(evi));
end
par.azm=azm;

for evi=1:length(azm)
    strazm(evi)={['bazm' num2str2(azm(evi),5,1)]};
end
set(handles.evlist,'string',strazm);

if par.g>length(par.st(par.h).evlo)
    par.g=length(par.st(par.h).evlo);
end

set(handles.evlist,'value',par.g);
par=plotdisp(par);
stmap=handles.stmap;
hold(stmap,'on');
PlotStationLocation(handles,par)
stlat=par.st(par.h).stla;
stlon=par.st(par.h).stlo;
plot(stmap,stlon,stlat,'ro');
hold(stmap,'off');

plotazm(par)

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
par.st(par.h)=[];
st=par.st;

if par.h>length(par.st)
    par.h=length(par.st);
else
    par.h=par.h;
end


for sti=1:length(st)
    stname(sti)=st(sti).stn;
    par.stloc(sti,:)=[st(sti).stla,st(sti).stlo];
    for i=1:length(par.periods);
        visom=st(sti).v0(:,i);
        visom(isnan(visom))=[];
        par.viso(sti,i)=median(visom);
    end
end

par.stname=stname;
set(handles.stlist,'string',stname);
set(handles.stlist,'value',par.h);

for i=1:length(par.periods);
    visom=par.viso(:,i);
    visom(isnan(visom))=[];
    par.vmall(i)=median(visom);
end


for evi=1:length(st(par.h).evla)
    azm(evi)=azimuth0(st(par.h).stla,st(par.h).stlo,st(par.h).evla(evi),st(par.h).evlo(evi));
end
par.azm=azm;

for evi=1:length(azm)
    strazm(evi)={['bazm' num2str2(azm(evi),5,1)]};
end
set(handles.evlist,'string',strazm);
set(handles.evlist,'value',1);


par=IMGstaDisp(par);
PlotStationLocation(handles,par)
stmp=handles.stmap;
hold(stmp,'on')
stlat=par.st(par.h).stla;
stlon=par.st(par.h).stlo;
plot(stmp,stlon,stlat,'ro');
hold(stmp,'off')

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
par.dvm=str2double(get(handles.paradvm,'string'));
par.dviso=str2double(get(handles.paradviso,'string'));
par.dvtance=str2double(get(handles.paradvtance,'string'));
par.dvazm=str2double(get(handles.paradvazm,'string'));
st=par.st;
save(par.outfile,'st','-v7.3')
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
sti=par.st;
v=sti(par.h).v0(par.g,:);
% v(ismember(par.periods,[par.PdlY par.Pvd par.Pvd2]))=nan;
v(ismember(par.periods,[par.PdlY]))=nan;

par.st(par.h).v0(par.g,:)=v;
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
sti=par.st;
v=sti(par.h).v0(par.g,:);
% v(ismember(par.periods,[par.PdlX par.Pvd par.Pvd2]))=nan;
v(ismember(par.periods,[par.PdlX]))=nan;

par.st(par.h).v0(par.g,:)=v;
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
plotazm(par);
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
stj=par.st(par.h);
if par.g<length(stj.evlo)
    par.g=par.g+1;
end
set(handles.evlist,'value',par.g);
par=plotdisp(par);
plotazm(par);
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
ds=par.ds;
hold(handles.figdsp1,'on')
par.dlX=[];
par.dlY=[];
v=ds.v;
periods=ds.period;

[y,x] = ginput(2);
deleteXi=find(v>min(x)&v<max(x));
deleteYi=find(periods>min(y)&periods<max(y));
scatter(periods(deleteXi),v(deleteXi),25,'fill','g')
hold(handles.figdsp1,'on')
scatter(periods(deleteYi),v(deleteYi),20,'fill','r')
hold(handles.figdsp1,'off')
par.PdlX=periods(deleteXi);
par.PdlY=periods(deleteYi);
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
par.dvm=str2double(get(handles.paradvm,'string'));
par.dviso=str2double(get(handles.paradviso,'string'));
par.dvtance=str2double(get(handles.paradvtance,'string'));
par.dvazm=str2double(get(handles.paradvazm,'string'));
loop=get(handles.Loops,'string');
loop=str2double(loop);
for loops=1:loop
    st=par.st;
    kh=0;
    iidg=[];
    for h=1:length(st)
        par.h=h;
        disp([num2str(h/length(st)*100),'%'])
        sti=st(h);
        for g=1:length(sti.evla) 
            par.g=g;
            vr=sti.v0(g,:);
            periods=par.periods;
            v=vr;
            prd=sti.periods(g,:);        
            idnan=~isnan(v);
            
            prd(isnan(v))=[];
            v(isnan(v))=[];
            
            idv1=abs(v-par.vmall(idnan));
            bad1= idv1>par.dvm;
            badpd1=prd(bad1);
            badid1=find(ismember(periods,badpd1));
            
            
            
            idv2=abs(v-par.viso(h,idnan));
            bad2= idv2>par.dviso;
            badpd2=prd(bad2);
            badid2=find(ismember(periods,badpd2));
                        



            if loops<3
                if length(v)>=25
                    pv = polyfit(prd,v,3);
                    vv=polyval(pv,prd);
                    iidg(h,g)=nan;
                elseif length(v)<25&&length(v)>10
                    pv = polyfit(prd,v,3);
                    vv=polyval(pv,prd);
                    iidg(h,g)=nan;
                else
                    iidg(h,g)=g;
                    continue
                end
            else
                pv = polyfit(prd,v,3);
                vv=polyval(pv,prd);
                iidg(h,g)=nan;
            end
                
            
            idv3=abs(v-vv);
            bad3= idv3>par.dvtance;
            badpd3=prd(bad3);
            badid3=find(ismember(periods,badpd3));
            
            %%
            azmLim=str2double(get(handles.azmLim,'string'));
            nbgrd=str2double(get(handles.nbgrd,'string'));
            hg=findnearDisp(par,nbgrd,azmLim);
            par.hg=hg;
            hi=hg.hi;
            gii=hg.gii;
            for i=1:length(hi)
                idh=hi(i);
                idg=cell2mat(gii(i));
                stgi=st(idh);
                vi(1:length(idg),:)=stgi.v0(idg,:);
%                 figure(999);
                for ip=1:length(vi(1,:))
                    vmip=vi(:,ip);
                    vmip(isnan(vmip))=[];
                    if ~isempty(vmip)
                        vmi(i,ip)=median(vmip);
                    else
                        vmi(i,ip)=nan;
                    end
                end
%                 plot(par.periods,vmi(i,:),'color',[0.4 0.4 0.4])
%                 hold on  
            end
          
            for ip=1:length(par.periods)
                vmip=vmi(:,ip);
                vmip(isnan(vmip))=[];
                vm(ip)=median(vmip);
            end
            prdvm=par.periods;
            prdvm(isnan(vm))=[];
            vm(isnan(vm))=[];
            
%             plot(prdvm,vm,'y','LineWidth',2)
%             idv4=abs(v-vm(ismember(prdvm,prd)));
%             bad4= idv4>par.dvazm;
%             badpd4=prd(bad4);
%             badid4=find(ismember(periods,badpd4));
            
            pvm = polyfit(prdvm,vm,3);
            vvm=polyval(pvm,prd);
%             plot(prd,vvm,'y','LineWidth',2)
            idv4=abs(v-vvm);
            bad4= idv4>par.dvazm;
            badpd4=prd(bad4);
            badid4=find(ismember(periods,badpd4));
            
            %%
%             plot(par.periods,par.viso(par.h,:),'g','LineWidth',2)
%             plot(par.periods,par.vmall,'r','LineWidth',2)
%             plot(par.periods,par.viso(par.h,:)-par.dvm,'k')
%             plot(par.periods,par.viso(par.h,:)+par.dvm,'k')
%             scatter(prd,v,80,'MarkerEdgeColor','k',...
%                 'MarkerFaceColor',[0 .7 .8])
            
            
%             scatter(badpd1,v(idv1>par.dvm),80,'MarkerEdgeColor','w',...
%                 'MarkerFaceColor','R')
%             scatter(badpd2,v(idv2>par.dviso),60,'MarkerEdgeColor','w',...
%                 'MarkerFaceColor','G')
%             scatter(badpd3,v(idv3>par.dvtance),40,'MarkerEdgeColor','w',...
%                 'MarkerFaceColor','k')
%             scatter(badpd4,v(idv4>par.dvazm),20,'MarkerEdgeColor','w',...
%                 'MarkerFaceColor','Y')
%             
%             xlim([min(par.periods)-5,max(par.periods)+5])
%             ylim([3,4.5])
%             drawnow;
%             hold off

            
            
            
            badid=[badid1 badid2 badid3 badid4];  
            badid=unique(badid);
            st(h).v0(g,badid)=nan;
            st(h).dv(g,badid)=nan;
            st(h).a0(g,badid)=nan;
            st(h).da(g,badid)=nan;
            st(h).r0(g,badid)=nan;
            st(h).dr(g,badid)=nan;
            st(h).g0(g,badid)=nan;
            st(h).dg(g,badid)=nan;
            st(h).vg(g,badid)=nan;
            st(h).azmo(g,badid)=nan;
            st(h).periods(g,badid)=nan;
            st(h).nsti(g,badid)=nan;
            

        end
    end
    
    
    for h=1:length(iidg(:,1))
        idg=iidg(h,:);
        idg(isnan(idg))=[];
        idg(idg==0)=[];
        st(h).evla(idg)=[];
        st(h).evlo(idg)=[];
        st(h).dpth(idg)=[];
        st(h).v0(idg,:)=[];
        st(h).dv(idg,:)=[];
        st(h).a0(idg,:)=[];
        st(h).da(idg,:)=[];
        st(h).r0(idg,:)=[];
        st(h).dr(idg,:)=[];
        st(h).g0(idg,:)=[];
        st(h).dg(idg,:)=[];
        st(h).vg(idg,:)=[];
        st(h).azmo(idg,:)=[];
        st(h).periods(idg,:)=[];
        st(h).nsti(idg,:)=[];
        if isempty(st(h).evla)
            kh=kh+1;
            idh(kh)=h;
        end
    end
    
    if kh>0
        st(idh)=[];
    end
end




par.st=st;
userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);
save(par.outfile,'st','-v7.3')
disp('done')







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
if par.h<length(par.st)
    par.h=par.h+1;
end
set(handles.stlist,'value',par.h);

stj=par.st(par.h);
for evi=1:length(stj.evla)
    azm(evi)=azimuth0(median(par.stloc(:,1)),median(par.stloc(:,2)),stj.evla(evi),stj.evlo(evi));
end
%    azm=sort(azm);
par.azm=azm;

for evi=1:length(azm)
    strazm(evi)={['azm' num2str2(azm(evi),5,1)]};
end
set(handles.evlist,'string',strazm);
set(handles.evlist,'value',1);
par=IMGstaDisp(par);
subplot(handles.stmap);
hold on
PlotStationLocation(handles,par)
stlat=par.st(par.h).stla;
stlon=par.st(par.h).stlo;
plot(stlon,stlat,'ro');
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
if par.h>1
    par.h=par.h-1;
end
set(handles.stlist,'value',par.h);

stj=par.st(par.h);
for evi=1:length(stj.evla)
    azm(evi)=azimuth0(median(par.stloc(:,1)),median(par.stloc(:,2)),stj.evla(evi),stj.evlo(evi));
end
azm=sort(azm);

for evi=1:length(azm)
    strazm(evi)={['azm' num2str2(azm(evi),5,1)]};
end
set(handles.evlist,'string',strazm);
set(handles.evlist,'value',1);
par=IMGstaDisp(par);
subplot(handles.stmap);
hold on
PlotStationLocation(handles,par)
stlat=par.st(par.h).stla;
stlon=par.st(par.h).stlo;
plot(stlon,stlat,'ro');
hold off
userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);

function PlotStationLocation(handles,par)
h1 = handles.stmap;
hold(h1, 'off');
g_allsta.lat=par.stloc(:,1);
g_allsta.lon=par.stloc(:,2);
plot(h1, g_allsta.lon, g_allsta.lat, 'k.');
hold(h1, 'on');
hg=par.hg;
hi=hg.hi;
plot(h1, g_allsta.lon(hi), g_allsta.lat(hi), 'r.');
plot(h1, g_allsta.lon(hg.h), g_allsta.lat(hg.h), 'b.')
set(h1,'Visible','off')



function par=IMGstaDisp(par)
handles=par.handles;
sti=par.st(par.h);
f1=handles.figdsp1;
hold(f1, 'off');

for i=1:length(sti.evla)
    periods=par.periods;
    v=sti.v0(i,:);
    periods(isnan(v))=[];
    v(isnan(v))=[];
    
    if ~isempty(v)
        
        if length(v)>25
            pv = polyfit(periods,v,3);
        elseif length(v)<=20&&length(v)>=10
            pv = polyfit(periods,v,3);
        else
            pv = polyfit(periods,v,2);
        end
        
%         vv=polyval(pv,periods);
        vv=v;
        azm=par.azm(i);
        
        if azm<0
            azm=360+azm;
        end
        
        if azm>180
            azm=azm-180;
        end
        
        
        if azm>=0&&azm<180/7
            plot(f1,periods,vv,'color',[255,0,0]/255)
        elseif azm>=180/7&&azm<2*180/7
            plot(f1,periods,vv,'color',[255,165,0]/255)
        elseif azm>=2*180/7&&azm<3*180/7
            plot(f1,periods,vv,'color',[255,255,0]/255)
        elseif azm>=3*180/7&&azm<4*180/7
            plot(f1,periods,vv,'color',[0,255,0]/255)
        elseif azm>=4*180/7&&azm<5*180/7
            plot(f1,periods,vv,'color',[0,255,255]/255)
        elseif azm>=5*180/7&&azm<6*180/7
            plot(f1,periods,vv,'color',[0,0,255]/255)
        else
            plot(f1,periods,vv,'color',[139,0,255]/255)
        end
        xlim(f1,[min(par.periods)-5,max(par.periods)+5])
        ylim(f1,[3,4.5])
        hold(f1, 'on');
    else
        plot(f1,[min(par.periods)-5,max(par.periods)+5],[4.5,4.5],'k')
    end
end


% id=flip(1:length(sti.evla));
% for ori=1:length(id)
%     set(handles.figdsp1.Children(ori),'tag',num2str((id(ori))))
% end
% hh=get(gca,'children');
% set(hh,'ButtonDownFcn','indxh=get(gcbo,''tag'');disp(str2double(indxh));yy=get(gcf,''userdata'');par=yy.par;handles=par.handles; set(handles.evlist,''value'',str2double(indxh)); par.handles=handles;userdata.par=par;set(gcf,''userdata'',userdata)')
% hold(f1, 'off');


function par=plotdisp(par)
handles=par.handles;
st=par.st;
periods=par.periods;
v=st(par.h).v0(par.g,:);
id=~isnan(v);
periods(isnan(v))=[];
v(isnan(v))=[];
dv2=abs(v-par.viso(par.h,id));
dv3=abs(v-par.vmall(id));

%%
azmLim=str2double(get(handles.azmLim,'string'));
nbgrd=str2double(get(handles.nbgrd,'string'));
hg=findnearDisp(par,nbgrd,azmLim);
par.hg=hg;
hi=hg.hi;
gii=hg.gii;
for i=1:length(hi)    
    idh=hi(i);
    idg=cell2mat(gii(i));
    sti=st(idh);
    vi(1:length(idg),:)=sti.v0(idg,:);
    f2=handles.figdsp1;
    for ip=1:length(vi(1,:))
        vmip=vi(:,ip);
        vmip(isnan(vmip))=[];
        if ~isempty(vmip)
            vmi(i,ip)=median(vmip);
        else
            vmi(i,ip)=nan;
        end
    end    
    plot(f2,par.periods,vmi(i,:),'g')
    hold(f2,'on')

end
%%
for ip=1:length(par.periods)
    vmipi=vmi(:,ip);    
    vmipi(isnan(vmipi))=[];
    vm(ip)=median(vmipi);
end
prdvm=par.periods;
prdvm(isnan(vm))=[];
vm(isnan(vm))=[];

pvm = polyfit(prdvm,vm,3);
vvm=polyval(pvm,periods);
plot(f2,periods,vvm,'b','LineWidth',3)
dv4=abs(v-vvm);   

f2=handles.figdsp1;
plot(f2,par.periods,par.viso(par.h,:),'r','LineWidth',3)
plot(f2,par.periods,par.vmall,'y','LineWidth',3)
plot(f2,par.periods,par.viso(par.h,:)-par.dvm,'b')
plot(f2,par.periods,par.viso(par.h,:)+par.dvm,'b')
scatter(f2,periods,v,40,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0 .7 .8])

if length(v)>25
    pv = polyfit(periods,v,3);
elseif length(v)<=25&&length(v)>=15
    pv = polyfit(periods,v,3);
else
    pv = polyfit(periods,v,3);
end
vv=polyval(pv,periods);
plot(f2,periods,vv,'k','LineWidth',1);
dv1=abs(v-vv);


scatter(f2,periods(dv1>par.dvtance),v(dv1>par.dvtance),20,'MarkerEdgeColor','w',...
    'MarkerFaceColor','R')
scatter(f2,periods(dv2>par.dviso),v(dv2>par.dviso),20,'MarkerEdgeColor','w',...
    'MarkerFaceColor','R')
scatter(f2,periods(dv3>par.dvm),v(dv3>par.dvm),20,'MarkerEdgeColor','w',...
    'MarkerFaceColor','R')
scatter(f2,periods(dv4>par.dvazm),v(dv4>par.dvazm),20,'MarkerEdgeColor','w',...
    'MarkerFaceColor','R')

xlim(f2,[min(par.periods)-5,max(par.periods)+5])
ylim(f2,[3,4.5])
hold(f2,'off')

ds.v=v;
ds.period=periods;
par.Pvd=periods(dv1>par.dvtance);
par.Pvd2=periods(dv2>par.dvm);
par.ds=ds;


function hg=findnearDisp(par,KmLim1,azmLim)
h=par.h;
g=par.g;
st=par.st;
stj=st(par.h);
for evi=1:length(stj.evla)
    azm(evi)=azimuth0(median(par.stloc(:,1)),median(par.stloc(:,2)),stj.evla(evi),stj.evlo(evi));
end
par.azm=azm;

azm0=azm(g);
if azm0<azmLim/2
    azmax=azm0+azmLim/2;
    azmin=360-azmLim/2+azm0;
    gi1=find(azm<=azmax&azm>=0);
    gi2=find(azm<=360&azm>=azmin);
    gi=[gi1,gi2];
elseif azm0>360-azmLim/2
    azmax=azmLim/2-360+azm0;
    azmin=azm0-azmLim/2;
    gi1=find(azm<=360&azm>=azmin);
    gi2=find(azm<=azmax>=0);
    gi=[gi1,gi2];
else
    azmax=azm0+azmLim/2;
    azmin=azm0-azmLim/2;
    gi=find(azm>=azmin&azm<=azmax);
end
hg.gi=gi;
hg.h=h;

la1=par.stloc(h,1);
lo1=par.stloc(h,2);
lamax=la1+KmLim1/deg2km(1);
lamin=la1-KmLim1/deg2km(1);
lomax=lo1+KmLim1/deg2km(1);
lomin=lo1-KmLim1/deg2km(1);
La=find(par.stloc(:,1)<=lamax&par.stloc(:,1)>=lamin);
Lo=find(par.stloc(:,2)<=lomax&par.stloc(:,2)>=lomin);
hi=intersect(La,Lo);
for i=1:length(hi)
    index=hi(i);
    sti=st(index);
    clear azmi
    for j=1:length(sti.evla)
        azmi(j)=azimuth0(sti.stla,sti.stlo,sti.evla(j),sti.evlo(j));
    end
    if azm0<azmLim/2
        azmax=azm0+azmLim/2;
        azmin=360-azmLim/2+azm0;
        gii1=find(azmi<=azmax&azmi>=0);
        gii2=find(azmi<=360&azmi>=azmin);
        gii(i)={[gii1,gii2]};
    elseif azm0>360-azmLim/2
        azmax=azmLim/2-360+azm0;
        azmin=azm0-azmLim/2;
        gii1=find(azmi<=360&azmi>=azmin);
        gii2=find(azmi<=azmax>=0);
        gii(i)={[gii1,gii2]};
    else
        azmax=azm0+azmLim/2;
        azmin=azm0-azmLim/2;
        gii(i)={find(azmi>=azmin&azmi<=azmax)};
    end
end

k=0;
for i=1:length(gii)
    flggi=cell2mat(gii(i));
    if isempty(flggi)
        k=k+1;
        idbad(k)=i;
    end  
end

if k>0
    hi(idbad)=[];
    gii(idbad)=[];
end
hg.hi=hi;
hg.gii=gii;



function plotazm(par)
handles=par.handles;
stj=par.st(par.h);
par.g=get(handles.evlist,'value');
par.h=get(handles.stlist,'value');
azmLim=str2double(get(handles.azmLim,'string'));
nbgrd=str2double(get(handles.nbgrd,'string'));
hg=findnearDisp(par,nbgrd,azmLim);
par.hg=hg;
evlo=stj.evlo(par.g);
evla=stj.evla(par.g);
stla=stj.stla;
stlo=stj.stlo;
af=handles.azmap;
%%world map
load coast
plot(af,long,lat,'Color',[150 150 150]/255)
hold(af,'on')
plot(af,stj.evlo,stj.evla,'.r')
scatter(af,evlo,evla,'b')
scatter(af,stlo,stla,40,'^','MarkerEdgeColor','k','MarkerFaceColor','y','LineWidth',1)
xlim(af,[-180 180])
ylim(af,[-90 90])

stgii=par.st(hg.hi);
gii=hg.gii;

for i=1:length(hg.hi)  
    plot(af,stgii(i).evlo(cell2mat(gii(i))),stgii(i).evla(cell2mat(gii(i))),'.g')    
end

stgi=par.st(hg.h);
gi=hg.gi;
plot(af,stgi.evlo(gi),stgi.evla(gi),'.g')

hold(af,'off')
set(af,'Visible','off')



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
