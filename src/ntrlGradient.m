%% [dvdx,dvdy]=ntrlGradient(lon,lat,value,dr,Method)
%        lon: x
%        lat: y
%        value: z
%        dr: dx=dy=dr
%        Method: must be one of the following:
%       'linear'  - (default) Linear interpolation
%       'nearest' - Nearest neighbor interpolation
%       'natural' - Natural neighbor interpolation 


function [dvdx,dvdy]=ntrlGradient(lon,lat,value,dr,Method)
y0=min(lat);
x0=min(lon);
dsx=ones(size(lon))*nan;
dsy=dsx;
for ist=1:length(lon)   
    dist=distance(y0,x0,lat(ist),lon(ist))*111.1949;
    azm=azimuth(y0,x0,lat(ist),lon(ist));
    dsx(ist)=sin(azm/180*pi)*dist;
    dsy(ist)=cos(azm/180*pi)*dist;
end
xyF = scatteredInterpolant(dsx,dsy,value,Method);
dsx1=dsx-111.1949*dr;
dsx2=dsx+111.1949*dr;
dsy1=dsy-111.1949*dr;
dsy2=dsy+111.1949*dr;
zx1=xyF(dsx1,dsy);
zx2=xyF(dsx2,dsy);
zy1=xyF(dsx,dsy1);
zy2=xyF(dsx,dsy2);
dvdx=(zx1+zx2-2*(value))/(2*111.1949*dr);
dvdy=(zy1+zy2-2*(value))/(2*111.1949*dr);
