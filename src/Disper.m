function varargout = Disper(varargin)
% DISPER MATLAB code for Disper.fig
%      DISPER, by itself, creates a new DISPER or raises the existing
%      singleton*.
%
%      H = DISPER returns the handle to a new DISPER or the handle to
%      the existing singleton*.
%
%      DISPER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISPER.M with the given input arguments.
%
%      DISPER('Property','Value',...) creates a new DISPER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Disper_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Disper_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Disper

% Last Modified by GUIDE v2.5 26-Aug-2021 18:45:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1; 
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Disper_OpeningFcn, ...
    'gui_OutputFcn',  @Disper_OutputFcn, ...
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


% --- Executes just before Disper is made visible.
function Disper_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Disper (see VARARGIN)

% Choose default command line output for Disper
handles.output = hObject;

% Update handles structure
% temp=load('./aniso/stationAniso/Aniso180allst.mat');
temp=load('./aniso/Aniso180allst.mat');
par.outfile='./aniso/Aniso180allst.mat';
par.dvm=0.6;
par.dv3=0.2;

Aniso=temp.Aniso;
period(1,:)=unique(Aniso.period);
Aniso.period=period;
for sti=1:length(Aniso.stn)
    stname(sti)={Aniso.stn(sti,:)};
end

par.handles=handles;
par.Aniso=Aniso;
par.stname=stname;
set(handles.dispList,'string',stname);
set(handles.dispList,'value',1);
par.h=get(handles.dispList,'value');
par=IMGstaDisp(par);
PlotStationLocation(handles,par)
userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);


function par=IMGstaDisp(par)
handles=par.handles;
Aniso=par.Aniso;
stname=par.stname;
for i=1:length(stname)
    periods=Aniso.period;
    v=Aniso.viso(i,:);
    if isfield(Aniso,'Qviso')
        vQ=Aniso.Qviso(i,:);
        periodsQ=Aniso.period;
    end
    
    periods(isnan(v))=[];
    v(isnan(v))=[];
    subplot(handles.Fig1);
    plot(periods,v,'.b');
    xlim([10,85])
    ylim([2.2,4.5])
    hold(handles.Fig1,'on')
    
    if isfield(Aniso,'Qviso')
        periodsQ(isnan(vQ))=[];
        vQ(isnan(vQ))=[];
        plot(periodsQ,vQ,'.r');
    end
    
end

for i=1:length(Aniso.period(1,:))
    if isfield(Aniso,'Qviso')
        vm=Aniso.Qviso(:,i);
    else
        vm=Aniso.viso(:,i);
    end
    vm(isnan(vm))=[];
    par.vm(i)=median(vm);
end

plot(Aniso.period,par.vm,'K')
plot(Aniso.period,par.vm-par.dvm,'g')
plot(Aniso.period,par.vm+par.dvm,'g')
hold(handles.Fig1,'off')




function par=plotdisp(par)
handles=par.handles;
Aniso=par.Aniso;
periods=Aniso.period;

if isfield(Aniso,'Qviso')
    v=Aniso.Qviso(par.h,:);
else
    v=Aniso.viso(par.h,:);
end

id=~isnan(v);
periods(isnan(v))=[];
v(isnan(v))=[];

if length(v)<5
    button=questdlg('points<5, remove or not','Note','yes','no','yes');
    if strcmp(button,'yes')
        Aniso.st(par.h,:)=[];
        Aniso.stn(par.h,:)=[];
        Aniso.viso(par.h,:)=[];
        Aniso.Fai(par.h,:)=[];
        Aniso.M(par.h,:)=[];
        Aniso.fitRsrm(par.h,:)=[];
        Aniso.a(par.h,:)=[];
        Aniso.b(par.h,:)=[];
        Aniso.stdv(par.h,:)=[];
        Aniso.stdv2(par.h,:)=[];
        Aniso.stdFai(par.h,:)=[];
        Aniso.stdM(par.h,:)=[];
        Aniso.stdA(par.h,:)=[];
        Aniso.stdB(par.h,:)=[];
        if isfield(Aniso,'Qviso')
            Aniso.Qviso(par.h,:)=[];
        end
        par.Aniso=Aniso;
        for sti=1:length(Aniso.stn)
            stname(sti)={Aniso.stn(sti,:)};
        end
        
        par.handles=handles;
        par.Aniso=Aniso;
        par.stname=stname;
        set(handles.dispList,'string',stname);
        set(handles.dispList,'value',par.h);
        par.h=get(handles.dispList,'value');
        PlotStationLocation(handles,par)
        return
    end
    
end


pv = polyfit(periods,v,3);
vv=polyval(pv,periods);
subplot(handles.Fig1);
scatter(periods,v,40,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0 .7 .8])
hold(handles.Fig1,'on')
xlim([10,85])
ylim([2.2,4.5])
plot(Aniso.period,par.vm,'r','LineWidth',1)
plot(periods,vv,'g','LineWidth',1);
plot(Aniso.period,par.vm-par.dvm,'b')
plot(Aniso.period,par.vm+par.dvm,'b')
dv=abs(v-vv);
dv2=abs(v-par.vm(id));
scatter(periods(dv>par.dv3),v(dv>par.dv3),20,'MarkerEdgeColor','w',...
    'MarkerFaceColor','R')
scatter(periods(dv2>par.dvm),v(dv2>par.dvm),20,'MarkerEdgeColor','w',...
    'MarkerFaceColor','R')

hold(handles.Fig1,'off')
subplot(handles.Map);
hold on
PlotStationLocation(handles,par)
stlat=Aniso.st(par.h,1);
stlon=Aniso.st(par.h,2);
plot(stlon,stlat,'ro');
hold off

ds.v=v;
ds.period=periods;
par.Pvd=periods(dv>par.dv3);
par.Pvd2=periods(dv2>par.dvm);
par.ds=ds;

function PlotStationLocation(handles,par)
h1 = handles.Map;
hold(h1, 'off');
Aniso=par.Aniso;
g_allsta.lat=Aniso.st(:,1);
g_allsta.lon=Aniso.st(:,2);
plot(h1, g_allsta.lon, g_allsta.lat, 'k.');
hold(h1, 'on');


% UIWAIT makes Disper wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Disper_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function dispList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dispList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in dispList.
function dispList_Callback(hObject, ~, handles) %#ok<DEFNU>
% hObject    handle to dispList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns dispList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dispList
userdata=get(gcf,'userdata');
par = userdata.par;
df0 = get(hObject,'String');
par.h=get(handles.dispList,'value');
par.h
disper0 = df0{get(hObject,'Value')};
par.disper0=disper0;
par=plotdisp(par);

userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);

% --- Executes on button press in Back.
function Back_Callback(hObject, eventdata, handles)
% hObject    handle to Back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par = userdata.par;
par.h=get(handles.dispList,'value');
if par.h>1
    par.h=par.h-1;
end
par=plotdisp(par);
set(handles.dispList,'value',par.h);
userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);


% --- Executes on button press in Next.
function Next_Callback(hObject, eventdata, handles)
% hObject    handle to Next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par = userdata.par;
par.h=get(handles.dispList,'value');
if par.h<length(par.stname)
    par.h=par.h+1;
end
par=plotdisp(par);
set(handles.dispList,'value',par.h);
userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);

% --- Executes on button press in Remove.
function Remove_Callback(hObject, eventdata, handles)
% hObject    handle to Remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par = userdata.par;
par.h=get(handles.dispList,'value');
Aniso=par.Aniso;

if isfield(Aniso,'Qviso')
    Aniso.Qviso(par.h,:)=[];
end

Aniso.viso(par.h,:)=[];
Aniso.a(par.h,:)=[];
Aniso.b(par.h,:)=[];
Aniso.M(par.h,:)=[];
Aniso.Fai(par.h,:)=[];
Aniso.stn(par.h,:)=[];
Aniso.st(par.h,:)=[];
Aniso.fitRsrm(par.h,:)=[];
Aniso.stdv(par.h,:)=[];
Aniso.stdA(par.h,:)=[];
Aniso.stdB(par.h,:)=[];
Aniso.stdM(par.h,:)=[];
Aniso.stdFai(par.h,:)=[];
par.Aniso=Aniso;
save(par.outfile,'Aniso')


for sti=1:length(Aniso.stn)
    stname(sti)={Aniso.stn(sti,:)};
end

if par.h>length(Aniso.st(:,1))
   par.h=length(Aniso.st(:,1));
end

par.stname=stname;
set(handles.dispList,'string',stname);
set(handles.dispList,'value',par.h);
PlotStationLocation(handles,par)
par=plotdisp(par);


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
par=plotdisp(par);

ds=par.ds;
hold(handles.Fig1,'on')
par.dlX=[];
par.dlY=[];
v=ds.v;
periods=ds.period;

[y,x] = ginput(2);
deleteXi=find(v>min(x)&v<max(x));
deleteYi=find(periods>min(y)&periods<max(y));
scatter(periods(deleteXi),v(deleteXi),25,'fill','g')
hold(handles.Fig1,'on')
scatter(periods(deleteYi),v(deleteYi),20,'fill','r')
hold(handles.Fig1,'off')
par.PdlX=periods(deleteXi);
par.PdlY=periods(deleteYi);
userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);

% --- Executes on button press in Dr.
function Dr_Callback(hObject, eventdata, handles)
% hObject    handle to Dr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par=userdata.par;
ds=par.ds;
Aniso=par.Aniso;

if isfield(Aniso,'Qviso')
    v=Aniso.Qviso(par.h,:);
else
    v=Aniso.viso(par.h,:);
end
% v(ismember(Aniso.period,[par.PdlY par.Pvd par.Pvd2]))=nan;
v(ismember(Aniso.period,[par.PdlY]))=nan;

Aniso.Qviso(par.h,:)=v;
par.Aniso=Aniso;
par=plotdisp(par);
save(par.outfile,'Aniso')

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
ds=par.ds;
Aniso=par.Aniso;

if isfield(Aniso,'Qviso')
    v=Aniso.Qviso(par.h,:);
else
    v=Aniso.viso(par.h,:);
end

% v(ismember(Aniso.period,[par.PdlX par.Pvd par.Pvd2]))=nan;
v(ismember(Aniso.period,[par.PdlX]))=nan;

Aniso.Qviso(par.h,:)=v;
par.Aniso=Aniso;
par=plotdisp(par);
save(par.outfile,'Aniso')

% Update handles structure
userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);

% --- Executes on button press in AotuPick.
function AotuPick_Callback(hObject, eventdata, handles)
% hObject    handle to AotuPick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(gcf,'userdata');
par = userdata.par;
Aniso=par.Aniso;
if isfield(Aniso,'Qviso')
    for i=1:length(Aniso.period);
        vmi=Aniso.Qviso(:,i);
        vmi(isnan(vmi))=[];
        vm(i)=median(vmi);
    end
    
    for i=1:length(Aniso.st(:,1))
        disp(i)
        clear vr id idp
        vr=Aniso.Qviso(i,:);
        prd=Aniso.period;
        id=~isnan(vr);
        nanid=isnan(vr);
        periods=prd;
        v=vr;
        periods(nanid)=[];
        v(nanid)=[];
        pv = polyfit(periods,v,3);
        vv=polyval(pv,periods);
        d1=abs(v-vv);
        d2=abs(v-vm(id));
        idv3=find(d1>par.dv3);
        idvm=find(d2>par.dvm);
        idv=[idv3 idvm];
        idv=unique(idv);
        
        if ~isempty(idv)
            for j=1:length(idv)
                idp(j)=find(Aniso.period==periods(idv(j)));
            end
        else
            idp=[];
        end
        
        vr(idp)=nan;
        Aniso.Qviso(i,:)=vr;
    end
    k=0;
    for ist=1:length(Aniso.st(:,1))
        vi=Aniso.Qviso(ist,:);
        vi(isnan(vi))=[];
        if length(vi)<5;
            k=k+1;
            badid(k)=ist;
        end
    end
    
    if k>0
        Aniso.Qviso(badid,:)=[];
        Aniso.viso(badid,:)=[];
        Aniso.a(badid,:)=[];
        Aniso.b(badid,:)=[];
        Aniso.M(badid,:)=[];
        Aniso.Fai(badid,:)=[];
        Aniso.stn(badid,:)=[];
        Aniso.st(badid,:)=[];
        Aniso.fitRsrm(badid,:)=[];
%         Aniso.stdv2(badid,:)=[];
        Aniso.stdv(badid,:)=[];
        Aniso.stdA(badid,:)=[];
        Aniso.stdB(badid,:)=[];
        Aniso.stdM(badid,:)=[];
        Aniso.stdFai(badid,:)=[];
    end
    
    save(par.outfile,'Aniso')
else
    for i=1:length(Aniso.period);
        vmi=Aniso.viso(:,i);
        vmi(isnan(vmi))=[];
        vm(i)=median(vmi);
    end
    
    for i=1:length(Aniso.st(:,1))
        clear vr id idp
        vr=Aniso.viso(i,:);
        prd=Aniso.period;
        
        id=~isnan(vr);
        nanid=isnan(vr);
        periods=prd;
        v=vr;
        
        periods(nanid)=[];
        v(nanid)=[];
        
        pv = polyfit(periods,v,4);
        vv=polyval(pv,periods);
        
        d1=abs(v-vv);
        d2=abs(v-vm(id));
        idv3=find(d1>par.dv3);
        idvm=find(d2>par.dvm);
        idv=[idv3 idvm];
        idv=unique(idv);
        
        if ~isempty(idv)
            for j=1:length(idv)
                idp(j)=find(Aniso.period==periods(idv(j)));
            end
        else
            idp=[];
        end
        vr(idp)=nan;
        Aniso.Qviso(i,:)=vr;
    end
    
    k=0;
    for ist=1:length(Aniso.st(:,1))
        vi=Aniso.Qviso(ist,:);
        vi(isnan(vi))=[];
        if length(vi)<5;
            k=k+1;
            badid(k)=ist;
        end
    end
    
    if k>0
        Aniso.Qviso(badid,:)=[];
        Aniso.viso(badid,:)=[];
        Aniso.a(badid,:)=[];
        Aniso.b(badid,:)=[];
        Aniso.M(badid,:)=[];
        Aniso.Fai(badid,:)=[];
        Aniso.stn(badid,:)=[];
        Aniso.st(badid,:)=[];
        Aniso.fitRsrm(badid,:)=[];
%         Aniso.stdv2(badid,:)=[];
        Aniso.stdv(badid,:)=[];
        Aniso.stdA(badid,:)=[];
        Aniso.stdB(badid,:)=[];
        Aniso.stdM(badid,:)=[];
        Aniso.stdFai(badid,:)=[];
    end
    
    
    save(par.outfile,'Aniso')
end

par.Aniso=Aniso;
for sti=1:length(Aniso.stn)
    stname(sti)={Aniso.stn(sti,:)};
end
par.stname=stname;
set(handles.dispList,'string',stname);
set(handles.dispList,'value',par.h);
% Update handles structure
userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);
