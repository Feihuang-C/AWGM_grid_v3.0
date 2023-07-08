function plot_map(wga,idx,cntrprd,evID)
suptitle(['ev=' evID '; T=' num2str(cntrprd) ' s'])
subplot(2,2,1)
dat=wga;
V=[dat.st(idx,1) dat.st(idx,2) dat.v(idx)];
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
title('£¨a£© velocity(km/s)')
hold off


subplot(2,2,2)
dat=wga;
V=[dat.st(idx,1) dat.st(idx,2) -dat.a0(idx)];
% VinGrd=gridatawg(lat,lon,V(:,1),V(:,2),V(:,3));
VinGrd = griddata(V(:,1),V(:,2),V(:,3),lat,lon);
[~,c]=contourf(lon,lat,VinGrd,200);
set(c,'edgecolor','none');
% c.LineColor='w';
colormap(flipud(jet))
caxis([-1*max([abs(min(VinGrd(:))),abs(max(VinGrd(:)))]),max([abs(min(VinGrd(:))),abs(max(VinGrd(:)))])]);
hold on

azm=azimuth(dat.ev(1,1),dat.ev(1,2),dat.st(:,1),dat.st(:,2));
azmo=azm-dat.a0;

a=sin(azm/180*pi);
b=cos(azm/180*pi);
ao=sin(azmo/180*pi);
bo=cos(azmo/180*pi);

ag=gridatawg(lat,lon,V(:,1),V(:,2),ao);
bg=gridatawg(lat,lon,V(:,1),V(:,2),bo);
agr=gridatawg(lat,lon,V(:,1),V(:,2),a);
bgr=gridatawg(lat,lon,V(:,1),V(:,2),b);

h1=quiver(lon,lat,ag,bg,'k');
h1=quiver(lon,lat,agr,bgr,'w');
shading interp
colorbar
box on
axis equal
ylim([min(V(:,1)),max(V(:,1))])
xlim([min(V(:,2)),max(V(:,2))])
title('(b) azimuth variation (deg)')
hold off


subplot(2,2,3)
dat=wga;
V=[dat.st(idx,1) dat.st(idx,2) dat.gs(idx)];
VinGrd = griddata(V(:,1),V(:,2),V(:,3),lat,lon);
% VinGrd=gridatawg(lat,lon,V(:,1),V(:,2),V(:,3));
% surface(lon,lat,VinGrd);
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
title('(c) geometrical spreadding (1/1000km)')
hold off


subplot(2,2,4)
dat=wga;
V=[dat.st(idx,1) dat.st(idx,2) dat.rd(idx)];
VinGrd = griddata(V(:,1),V(:,2),V(:,3),lat,lon);
% VinGrd=gridatawg(lat,lon,V(:,1),V(:,2),V(:,3));

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
title('(d) radiation pattern(1/deg)')
hold off




function zg=gridatawg(xg,yg,x,y,z)
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