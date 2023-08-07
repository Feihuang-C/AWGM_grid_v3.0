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

% Last Modified by GUIDE v2.5 10-Jun-2023 16:01:20

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
function Disper_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Disper (see VARARGIN)

% Choose default command line output for Disper
handles.output = hObject;

% Update handles structure
% temp=load('./aniso/stationAniso/Aniso180allst.mat');
[file,path,indx] = uigetfile('*.mat');
if indx==0;
    return
end
infile=[path file];
temp=load(infile);
Aniso=temp.Aniso;
strinf=strsplit(file,'_');
if strcmp(cell2mat(strinf(1)),'AApick')
   par.outfile=[path file];
else
   par.outfile=[path 'AApick_' file];
end
par.dvm=str2double(get(handles.vlimt1,'string'));
par.dv3=str2double(get(handles.vlimt2,'string'));
save(par.outfile,'Aniso')
period(1,:)=unique(Aniso.period);
Aniso.period=period;
stname(1:length(Aniso.st(:,1)))={nan};
for sti=1:length(Aniso.st(:,1))
    stname(sti)={num2str(sti)};
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
 
    periods(isnan(v))=[];
    v(isnan(v))=[];
    subplot(handles.Fig1);
    plot(periods,v,'.b');
    xlim([min(periods)-5,max(periods)+5])
    ylim([2,5])
    hold(handles.Fig1,'on')    
    
end
xlabel('T (s)')
ylabel('v (km/s)')
for i=1:length(Aniso.period(1,:))
    vm=Aniso.viso(:,i);
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
prd=Aniso.period;
periods=prd;
viso=Aniso.viso(par.h,:);
v=viso;
id=~isnan(v) & v<10;
periods(isnan(v) | v>=10)=[];
v(isnan(v) | v>=10)=[];
pv = polyfit(periods,v,3);
vv=polyval(pv,prd);
subplot(handles.Fig1);
scatter(periods,v,40,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0 .7 .8])
hold(handles.Fig1,'on')
xlabel('T (s)')
ylabel('v (km/s)')
box on
xlim([min(periods)-5,max(periods)+5])
ylim([2,5])
plot(Aniso.period,par.vm,'r','LineWidth',1)
plot(prd,vv,'g','LineWidth',1);
plot(Aniso.period,par.vm-par.dvm,'b')
plot(Aniso.period,par.vm+par.dvm,'b')
dv=abs(viso-vv);
dv2=abs(viso-par.vm);
scatter(prd(dv>par.dv3),viso(dv>par.dv3),20,'MarkerEdgeColor','w',...
    'MarkerFaceColor','R')
scatter(prd(dv2>par.dvm),viso(dv2>par.dvm),20,'MarkerEdgeColor','w',...
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
par.Pvd=prd(dv>par.dv3);
par.Pvd2=prd(dv2>par.dvm);
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
function dispList_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
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

stname(1:length(Aniso.st(:,1)))={' '};
for sti=1:length(Aniso.st(:,1))
    stname(sti)={num2str(sti)};
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
v=Aniso.viso(par.h,:);
% v(ismember(Aniso.period,[par.PdlY par.Pvd par.Pvd2]))=nan;
Aniso.vori(par.h,:)=v;
v(ismember(Aniso.period,[par.PdlY]))=nan;

Aniso.viso(par.h,:)=v;
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


v=Aniso.viso(par.h,:);
v(ismember(Aniso.period,[par.PdlX]))=nan;
Aniso.viso(par.h,:)=v;
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

vm(1:length(Aniso.period))=nan;
for i=1:length(Aniso.period);
    vmi=Aniso.viso(:,i);
    vmi(isnan(vmi))=[];
    vm(i)=median(vmi);
end

k=0;
for i=1:length(Aniso.st(:,1))
    disp(i)
    clear viso id idp
    viso=Aniso.viso(i,:);
    prd=Aniso.period;
    
    id=~isnan(viso);
    nanid=isnan(viso);
    periods=prd;
    v=viso;
    
    periods(nanid)=[];
    v(nanid)=[];
    
    pv = polyfit(periods,v,3);
    d1=abs(viso-polyval(pv,prd));
    d2=abs(viso-vm);
    idv3=d1>par.dv3;
    idvm=d2>par.dvm;
    idp=idv3 | idvm;   
    viso(idp)=nan;  
    if length(viso(~isnan(viso)))<5   
        k=k+1;
        badid(k)=i;
        Aniso.viso(i,:)=nan;
    else    
        Aniso.viso(i,:)=viso;
    end
end
if k>0
Aniso.st(badid,:)=[];
Aniso.viso(badid,:)=[];
Aniso.a(badid,:)=[];
Aniso.b(badid,:)=[];
Aniso.M(badid,:)=[];
Aniso.Fai(badid,:)=[];
Aniso.st(badid,:)=[];
Aniso.fitRsrm(badid,:)=[];
Aniso.stdv(badid,:)=[];
Aniso.stdA(badid,:)=[];
Aniso.stdB(badid,:)=[];
Aniso.stdM(badid,:)=[];
Aniso.stdFai(badid,:)=[];
end

save(par.outfile,'Aniso')


par.Aniso=Aniso;
stname(1:length(Aniso.st))={''};
for sti=1:length(Aniso.st(:,1))
    stname(sti)={num2str(sti)};
end
par.stname=stname;
set(handles.dispList,'string',stname);
set(handles.dispList,'value',par.h);
% Update handles structure
userdata.par=par;
set(gcf,'userdata',userdata);
guidata(hObject, handles);



function vlimt1_Callback(hObject, eventdata, handles)
% hObject    handle to vlimt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vlimt1 as text
%        str2double(get(hObject,'String')) returns contents of vlimt1 as a double
userdata=get(gcf,'userdata');
par = userdata.par;
par.dvm=str2double(get(handles.vlimt1,'string'));
userdata.par=par;
set(gcf,'userdata',userdata);

% --- Executes during object creation, after setting all properties.
function vlimt1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vlimt1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vlimt2_Callback(hObject, eventdata, handles)
% hObject    handle to vlimt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vlimt2 as text
%        str2double(get(hObject,'String')) returns contents of vlimt2 as a double
userdata=get(gcf,'userdata');
par = userdata.par;
par.dv3=str2double(get(handles.vlimt2,'string'));
userdata.par=par;
set(gcf,'userdata',userdata);

% --- Executes during object creation, after setting all properties.
function vlimt2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vlimt2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
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
    case 'f'
        Next_Callback(handles.Next,[],handles)
    case 'd'
        Back_Callback(handles.Back,[],handles)
    case 'r'
        Dr_Callback(handles.Dr,[],handles)
    case 'g'
        Dg_Callback(handles.Dg,[],handles)
    otherwise      
end

function Fig1_ButtonDownFcn(hObject, eventdata, handles)

% switch (get(gcbf,'SelectionType'))
%     case 'alt'
%         Pick_Callback(handles.Pick,[],handles)
%     otherwise
% end

function figure1_WindowButtonDownFcn(hObject, eventdata, handles)

switch (get(gcbf,'SelectionType'))
    case 'alt'
        Pick_Callback(handles.Pick,[],handles)
    case 'normal'
%         Next_Callback(handles.Next,[],handles)
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
    Back_Callback(handles.Next,[],handles)    
elseif count==1
    Next_Callback(handles.Next,[],handles)    
end




function figure1_SizeChangedFcn(hObject, eventdata, handles)