function sks_pms_ani
close all
clear rall

pthi0=['..\..\data\sws\'];
flpms0=[pthi0 'ss_pms_all.txt'];
flpms1=[pthi0 'ss_pms_chenyun_etibet.txt']
flpms2=[pthi0 'ss_pms_sunya-rf.txt']

pthi1=['..\..\data\gps\'];
flgps0=[pthi1 'GPS_vel_all_EuraRef_vect.txt'];
flgps1=[pthi1 'GPS_shen_zhen_kang_vect.txt']

%plot_pms_vs_pms_statistics(flpms1,flpms2,200)

%%compare PMS .vs. NZ
flnz=['..\..\run03\t01i02\vslay_dep0000_0050_da10_ans3d_smth.dat'];
%flnz=['..\run01\inv03\vslay_dep0000_0050_da10_ans3d_smth.dat'];

ddir0=plot_pms_vs_nz_statistics(flnz,flpms0, 300);
%ddir1=plot_pms_vs_nz_statistics(flnz,flpms1, 400);
%ddir2=plot_pms_vs_nz_statistics(flnz,flpms2, 500);

ddir2=plot_pms_vs_nz_statistics(flnz,flgps0, 600);
ddir2=plot_pms_vs_nz_statistics(flnz,flgps1, 700);

return

ddir=plot_pms_vs_gps_statistics(flgps0,flpms1,800);
ddir=plot_pms_vs_gps_statistics(flgps0,flpms2,900);


%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%
function ddir=plot_pms_vs_nz_statistics(flnz,flpms,ifig)
%%plot dsp variation with azimuths/prd in one block
[stn,lat,lon,fpo,aoa,f1,f2,a1,a2,nst]=textread(flpms,'%s %f %f %f %f %f %f %f %f %f');
stn2=unique(stn);
ns2=length(stn2);
for ii=1:ns2
    idx=find(strcmp(stn2(ii),stn));
    lat1(ii)=lat(idx(1));
    lon1(ii)=lon(idx(1));
    fp0=median(fpo(idx));
    if(fp0<0)
        fp0=fp0+180;
    elseif(fp0>180)
        fp0=fp0-180;
    end
    fpo1(ii)=fp0;    
end

[lat2,lon2,vmean,aamp,fpo2,kbk]=textread(flnz,'%f %f %f %f %f %f');

na1=length(lat1);
na2=length(lat2);
ncom=0;
dmax=0.5;
for ib=1:na2 
    xi=lon2(ib);
    yi=lat2(ib);
    fp0=fpo2(ib);
    if(fp0<0)
        fp0=fp0+180;
    elseif(fp0>180)
        fp0=fp0-180;
    end
    dx=abs(lon1-xi);
    dy=abs(lat1-yi);
    idx1=(dx<dmax & dy<dmax);
    idx=find(idx1);
    
    %[ddi,i]=min(dd(idx));
    %idx=idx(i);
    if(~isempty(idx))
        fp1=fpo1(idx);
        dis=sqrt(dx(idx).*dx(idx)+dy(idx).*dy(idx));
        [dis0,ii]=min(dis);
        dd=abs(fp1(ii)-fp0);   %%values from the closest grids
        %dd=min(abs(fp1-fp0));  %%minimum value within 0.5 degree grids
        if(dd>90)
            dd=180-dd;
        end
        ncom=ncom+1;
        ddir(ncom)=dd;
    end
end
figure(ifig)
%[nc xc]=hist(ddir,[0:10:180]);
subplot(1,2,1)
rose(ddir/180*pi,24);
view(90,-90)
subplot(1,2,2)
rose(ddir/180*pi,36);
%title(['Med=' num2str(median(ddir))])
view(90,-90)
msg=[ncom median(ddir) mean(ddir)]



%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%
function ddir=plot_pms_vs_gps_statistics(flgps,flpms,ifig)
%%plot dsp variation with azimuths/prd in one block
[stn,lat,lon,fpo,aoa,f1,f2,a1,a2,nst]=textread(flpms,'%s %f %f %f %f %f %f %f %f %f');
stn2=unique(stn);
ns2=length(stn2);
for ii=1:ns2
    idx=find(strcmp(stn2(ii),stn));
    lat1(ii)=lat(idx(1));
    lon1(ii)=lon(idx(1));
    fp0=median(fpo(idx));
    if(fp0<0)
        fp0=fp0+180;
    elseif(fp0>180)
        fp0=fp0-180;
    end
    fpo1(ii)=fp0;    
end

[stn2,lat2,lon2,fpo2,aamp,kbk,a,b]=textread(flgps,'%s %f %f %f %f %f %f %f');

na1=length(lat1);
na2=length(lat2);
ncom=0;
dmax=0.25;
for ib=1:na2 
    xi=lon2(ib);
    yi=lat2(ib);
    fp0=fpo2(ib);
    dx=abs(lon1-xi);
    dy=abs(lat1-yi);
    idx1=(dx<dmax & dy<dmax);
    idx=find(idx1);
    
    %[ddi,i]=min(dd(idx));
    %idx=idx(i);
    if(~isempty(idx))
        fp1=fpo1(idx);
        dis=sqrt(dx(idx).*dx(idx)+dy(idx).*dy(idx))
        ii=min(dis);
        dd=abs(fp1(ii)-fp0);
        %dd=min(abs(fp1-fp0));
        if(dd>90)
            dd=180-dd;
        end
        ncom=ncom+1;
        ddir(ncom)=dd;
    end
end
figure(ifig)
%[nc xc]=hist(ddir,[0:10:180]);
subplot(1,2,1)
rose(ddir/180*pi,24);
view(90,-90)
subplot(1,2,2)
rose(ddir/180*pi,36);
%title(['Med=' num2str(median(ddir))])
view(90,-90)
msg=[ncom median(ddir) mean(ddir)]



%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%
function plot_pms_vs_pms_statistics(flpms1,flpms2,ifig)
%%plot dsp variation with azimuths/prd in one block
[stn1,lat1,lon1,fpo1,aoa1,a,b,c,d,f]=textread(flpms1,'%s %f %f %f %f %f %f %f %f %f');
stn1u=unique(stn1);
ns1u=length(stn1u)
for ii=1:ns1u
    idx=find(strcmp(stn1u(ii),stn1));
    lat1u(ii)=lat1(idx(1));
    lon1u(ii)=lon1(idx(1));
    fp0=median(fpo1(idx));
    if(fp0<0)
        fp0=fp0+180;
    elseif(fp0>180)
        fp0=fp0-180;
    end
    fpo1u(ii)=fp0;    
end

[stn2,lat2,lon2,fpo2,aoa2]=textread(flpms2,'%s %f %f %f %f');
stn2u=unique(stn2);
ns2u=length(stn2u)
for ii=1:ns2u
    idx=find(strcmp(stn2u(ii),stn2));
    lat2u(ii)=lat2(idx(1));
    lon2u(ii)=lon2(idx(1));
    fp0=median(fpo2(idx));
    if(fp0<0)
        fp0=fp0+180;
    elseif(fp0>180)
        fp0=fp0-180;
    end
    fpo2u(ii)=fp0;    
end

ncom=0;
for ib=1:ns2u 
    xi=lon2u(ib);
    yi=lat2u(ib);
    fp0=fpo2u(ib);
    dx=lon1u-xi;
    dy=lat1u-yi;
    dd=sqrt(dx.*dx+dy.*dy);
    stni=stn2u(ib);
    idx=find(strcmp(stni, stn1u));
    if(~isempty(idx))
        fp1=median(fpo1u(idx));
        ncom=ncom+1;
        dd=min(abs(fpo1u(idx)-fp0));
        if(dd>90)
            dd=180-dd;
        end
        
        ddir(ncom)=dd;
        
    end       
end
figure(ifig);
subplot(1,2,1)
%[nc xc]=hist(ddir,[0:10:180]);
rose(ddir/180*pi,24);
view(90,-90);
msg1=[ncom median(ddir) mean(ddir)]



%%
%%
%%
function sks_ani(fsks,fout)
[stn,lat,lon,fpo,aoa]=textread(fsks,'%s\t%f\t%f\t%f\t%f\n');
stn2=unique(stn);

ns=length(stn)
ns2=length(unique(stn))

fid=fopen(fout,'w')
for is=1:ns2
    idx=find(strcmp(stn2(is),stn));    
    lat2(is)=lat(idx(1));
    lon2(is)=lon(idx(1));
    fpo2(is)=mean(fpo(idx));
    aoa2(is)=mean(aoa(idx));
    
    fprintf(fid,'%s %f %f %f %f\n',cell2mat(stn2(is)),lat2(is),lon2(is),fpo2(is),aoa2(is));
end
fclose(fid)

%%
[nc xc]=hist(aoa2,[0:0.2:2]);
bar(xc,nc);


%%
%% pms
%%
function pms_ani(fli,ifig)
[stn,lat,lon,fpo,aoa,a,b,c,d,f]=textread(fli,'%s %f %f %f %f %f %f %f %f %f');
stn2=unique(stn);

ns=length(stn)
ns2=length(unique(stn))
%%
figure(ifig)
[nc xc]=hist(aoa,[0:0.1:1.6]);
bar(xc,nc);
xlabel('relative time ')





