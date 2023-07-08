function plot_map2(wga,cntrprd,evID)
suptitle(['ev=' evID '; T=' num2str(cntrprd) ' s'])
subplot(2,2,1)
dat=wga;
V=[dat.st(:,1) dat.st(:,2) dat.v];
lat = unique((V(:,1)));
lon = unique((V(:,2)));
[lat,lon] = meshgrid(lat,lon);
% VinGrd=gridatawg(lat,lon,V(:,1),V(:,2),V(:,3));
VinGrd = griddata(V(:,1),V(:,2),V(:,3),lat,lon);
[~,c]=contourf(lon,lat,VinGrd,200);
set(c,'edgecolor','none');
% c.LineColor='w';
colormap(flipud(jet))
hold on
% caxis([-4000 4000]);
shading interp
colorbar
box on
axis equal
ylim([min(V(:,1)),max(V(:,1))])
xlim([min(V(:,2)),max(V(:,2))])
title('(a) v_s_m_r_v (km/s)')
hold off

subplot(2,2,2)
dat=wga;
V=[dat.st(:,1) dat.st(:,2) dat.va];
% VinGrd=gridatawg(lat,lon,V(:,1),V(:,2),V(:,3));
VinGrd = griddata(V(:,1),V(:,2),V(:,3),lat,lon);
[~,c]=contourf(lon,lat,VinGrd,200);
set(c,'edgecolor','none');
% c.LineColor='w';
colormap(flipud(jet))
hold on
% caxis([-4000 4000]);
shading interp
colorbar
box on
axis equal
ylim([min(V(:,1)),max(V(:,1))])
xlim([min(V(:,2)),max(V(:,2))])
title('(b) v_s_m_r_a (km/s)')
hold off


subplot(2,2,3)
dat=wga;
V=[dat.st(:,1) dat.st(:,2) dat.vcr];
% VinGrd=gridatawg(lat,lon,V(:,1),V(:,2),V(:,3));

VinGrd = griddata(V(:,1),V(:,2),V(:,3),lat,lon);
[~,c]=contourf(lon,lat,VinGrd,200);
set(c,'edgecolor','none');
% c.LineColor='w';
colormap(flipud(jet))
hold on
% caxis([-4000 4000]);
shading interp
colorbar
box on
axis equal
ylim([min(V(:,1)),max(V(:,1))])
xlim([min(V(:,2)),max(V(:,2))])
title('(c) vcr_s_m_r_v (km/s)')
hold off

subplot(2,2,4)
dat=wga;
V=[dat.st(:,1) dat.st(:,2) dat.dvcr];
% VinGrd=gridatawg(lat,lon,V(:,1),V(:,2),V(:,3));

VinGrd = griddata(V(:,1),V(:,2),V(:,3),lat,lon);
[~,c]=contourf(lon,lat,VinGrd,200);
set(c,'edgecolor','none');
% c.LineColor='w';
colormap(flipud(jet))
hold on
% caxis([-4000 4000]);
shading interp
colorbar

azm=azimuth(dat.evla,dat.evlo,dat.st(:,1),dat.st(:,2));
azmo=azm-dat.az;
a=sin(azmo/180*pi);
b=cos(azmo/180*pi);
ag=gridatawg(lat,lon,V(:,1),V(:,2),a);
bg=gridatawg(lat,lon,V(:,1),V(:,2),b);
h1=quiver(lon,lat,ag,bg,'k');
shading interp
colorbar
box on
axis equal
ylim([min(V(:,1)),max(V(:,1))])
xlim([min(V(:,2)),max(V(:,2))])
title('(d) dvcr (km/s)')
hold off

function zg=gridatawg(yg,xg,y,x,z)
[n,m]=size(xg);
for isg=1:n
    for jsg=1:m
        idx=find(x==xg(isg,jsg) & y==yg(isg,jsg));
        if isempty(idx)
            zg(isg,jsg)=nan;
        else
            zg(isg,jsg)=z(idx);
        end
        
    end
end
