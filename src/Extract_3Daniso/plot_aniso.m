dat=vmod;
x=[dat.xbk];
y=[dat.ybk];
xg=unique(x);
yg=unique(y);
yalim=[26 32];
yg(yg>yalim(2)|yg<yalim(1))=[];
[xg,yg] = meshgrid(xg,yg);
fi=0;
DEPTH=dat(1).dp;
DEP=[10 20 40 60 90];
visos=cell2mat({dat.viso}');
FPDs=cell2mat({dat.FPD}');
MOAs=cell2mat({dat.MOA}');
for i=DEP
    fi=fi+1;
    ids=find(DEPTH==i);
    Visos=visos(:,ids);
    FPD=FPDs(:,ids);
    MOA=MOAs(:,ids);
    MOA(MOA<0)=0;
    MOA(MOA>5)=5;
    a=MOA.*sin(FPD/180*pi);
    b=MOA.*cos(FPD/180*pi);
    [n,m]=size(xg);
    MOA(Visos<3)=0;
    Visos(Visos<3)=nan;
%     for isg=1:n
%         for jsg=1:m
%             idx=find(x==xg(isg,jsg) & y==yg(isg,jsg));
%             if isempty(idx)
%                 vGrd(isg,jsg)=nan;
%                 ag(isg,jsg)=nan;
%                 bg(isg,jsg)=nan;
%             else
%                 vGrd(isg,jsg)=Visos(idx);
%                 ag(isg,jsg)=a(idx);
%                 bg(isg,jsg)=b(idx);
%             end
%         end
%     end
            vGrd=griddata(x,y,Visos,xg,yg);
            ag=griddata(x,y,a,xg,yg);
            bg=griddata(x,y,b,xg,yg);
    
    figure
    f=contourf(xg,yg,vGrd);
    colormap(flipud(jet))
    hold on
    h1=quiver(xg,yg,ag,bg,'k');
    h1.ShowArrowHead = 'off';
    h2=quiver(xg,yg,-ag,-bg,'k');
    h2.ShowArrowHead = 'off';
    hold off
    axis equal
end