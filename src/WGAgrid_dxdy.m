%% ---comput spatial gradient for a subarray----
%--------------------------------------------
function [dudx,dudy,stc,uprdt]=WGAgrid_dxdy(st,ev,fband,par,~)
vavg=par.vavg;
damp=0.01;
fc=mean(fband);
azpth=azimuth0(st(1).st(1),st(1).st(2),ev(1),ev(2))+180;
azpth(azpth<0)=azpth(azpth<0)+360;
azpth(azpth>360)=azpth(azpth>360)-360;
stc=st(1);
nst=length(st);
[uobs,dx,dy,wi,G,~]=geowiF(st,nst,azpth,2,vavg,damp,fc,par);
loopflg=1;
loop=0;
while loopflg==1 && loop<2
    loop=1+loop;
    W=diag(wi,0);
    if par.kwi>0
        u0=pinv(G'*W*G)*G'*W*uobs;
    else
        u0=pinv(G)*uobs;
    end
    stc(1).dat=u0(1,:);
    dudx=u0(2,:);
    dudy=u0(3,:);
    uprdt=ones(size(uobs));
    Rst=ones(1,nst-1);
    g=0;
    for ist = 1:nst-1
        uprdt(ist,:)=u0(1,:)+dudx*dx(ist)+dudy*dy(ist);
        Rij = corrcoef(uprdt(ist,:), uobs(ist,:));
        Rst(ist)=Rij(1,2);
        if Rst(ist)<par.Rcof
            wi(ist)=wi(ist)*Rst(ist)^2;
            g=g+1;            
        end
    end
    if g>0
        loopflg=1;
    else
        loopflg=0;        
    end
end
 
return
%% return
figure(956) %#ok<UNRCH>
U=zeros(1,length(uobs(1,:)));
for ist = 1:length(dx)
    uprdt(ist,:)=u0(1,:)+dudx*dx(ist)+dudy*dy(ist);
    U=U+abs(uprdt(ist,:));
    plot(uobs(ist,:)+max(U)*3+2,'b-') 
    hold on
    plot(uprdt(ist,:)+max(U)*3+2,'--r')
    text(5,max(U)*3+0.5*max(uobs(ist,:)),num2str(Rst(ist)));
end
hold off



function [ui,dx,dy,wi,G,er]=geowiF(st,nst,azpth,sc1,vavg,damp,fc,par)
np=length(st(sc1).dat);
ui=ones(nst-1,np)*nan;
dx=ones(nst-1,1)*nan;
dy=ones(nst-1,1)*nan;
wi=ones(nst-1,1)*nan;
er=ones(nst-1,1)*nan;
for is=2:nst
    ui(is-1,1:np)=st(is).dat;    
    dxi=(st(is).st(2)-st(1).st(2));
    dyi=(st(is).st(1)-st(1).st(1));
    if dxi==0
        dx(is-1)=0;
    else
        dx(is-1)=distance0(st(1).st(1),st(1).st(2),st(1).st(1),st(is).st(2))*111.1949*dxi/abs(dxi);
    end
    if dyi==0
        dy(is-1)=0; 
    else
        dy(is-1)=distance0(st(1).st(1),st(1).st(2),st(is).st(1),st(1).st(2))*111.1949*dyi/abs(dyi);
    end
    wi(is-1)=compweight(azpth,st(1),st(is),damp,par.kwi);
    er(is-1)=CompErrors(azpth,st(1),st(is),vavg,fc);    
end
G=ones(nst-1,3);
G(:,2)=dx;
G(:,3)=dy;






%%    
function [kmax]=findmaxima(dat,mmax) %#ok<DEFNU>
diff0=diff(dat);
ndf=length(diff0);
nmax=0;
pmax=ones((ndf-1),1)*nan;
for id=2:ndf-1
    if diff0(id)==0 && diff0(id-1)>0 && diff0(id+1)<0
        nmax=nmax+1;
        pmax(nmax)=id;
        %[blw(nmax),bup(nmax)]=errbars(dat,id,0.95);
    elseif diff0(id)>0 && diff0(id+1)<0
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
kmax=ones(mmax,1)*nan;
for ii=1:mmax
    [~,kk]=max(dat(pmax));
    kmax(ii)=pmax(kk);
    pmax(kk)=[];
end

%%
%% Calculate the weights for station pairs.
%%
function wi=compweight(azpth,st1,st2,damp,kwi)
if ~isempty(azpth)
    dst=distance0(st1.st(1),st1.st(2),st2.st(1),st2.st(2))*111.1949;
    if dst==0
        wi=0.1/damp;
        wi=wi^kwi;
    else
        azst=azimuth0(st1.st(1),st1.st(2),st2.st(1),st2.st(2));
        daz=(azpth-azst)/180*pi;
        wi=dst*cos(daz);
        wi=0.1/(abs(wi)+damp);
        wi=wi^kwi;
    end
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
    lamda=vel/fc/111.1949;
    delta=pi/lamda*abs(wi);
else
    delta=[];
end

