%---------azimuth between two points-------------
function azm=azimuth0(ela,elo,sla,slo)
re=1.0;
pp=pi/180.0;
dis=distance0(ela,elo,sla,slo)*pp;
ela=ela*pp;elo=elo*pp;sla=sla*pp;slo=slo*pp;
xe=re*[cos(ela)*cos(elo) cos(ela)*sin(elo) sin(ela)];
xs=re*[cos(sla)*cos(slo) cos(sla)*sin(slo) sin(sla)];
b=cross(xs,xe)/re/re/sin(dis);
ula=[cos(elo)*sin(ela) sin(elo)*sin(ela) -cos(ela)];
ulo=[-sin(elo) cos(elo) 0];
cosazm=b*ulo';
sinazm=b*ula';
azm=atan2(sinazm,cosazm)/pp;
if azm<0.0
    azm=azm+360; 
end

function dis=distance0(ela,elo,sla,slo)
re=1.0;
pp=pi/180.0;
ela=ela*pp;elo=elo*pp;sla=sla*pp;slo=slo*pp;
xe=re*[cos(ela)*cos(elo) cos(ela)*sin(elo) sin(ela)];
xs=re*[cos(sla)*cos(slo) cos(sla)*sin(slo) sin(sla)];
cosdis=(xe*xs')/re/re;
dis=acos(cosdis);
dis=dis/pp;

