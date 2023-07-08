%---------distance between two points--------------
%-ela,sla are latitudes of earthquake and stations-
%-elo,slo are longitudes of earthquake and station-
%-------all angle inputs are in degrees------------
function dis=distance0(ela,elo,sla,slo)
re=1.0;
pp=pi/180.0;
ela=ela*pp;elo=elo*pp;sla=sla*pp;slo=slo*pp;
xe=re*[cos(ela)*cos(elo) cos(ela)*sin(elo) sin(ela)];
xs=re*[cos(sla)*cos(slo) cos(sla)*sin(slo) sin(sla)];
cosdis=(xe*xs')/re/re;
dis=acos(cosdis);
dis=dis/pp;

