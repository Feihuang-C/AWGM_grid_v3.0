    Aniso.period=10:120;
    x=Aniso.st(:,2);
    y=Aniso.st(:,1);
    load('gstall.mat')
    xg=unique(x); 
    yg=unique(y);
    yalim=[21 45];
    yg(yg>yalim(2)|yg<yalim(1))=[];
    [xg,yg] = meshgrid(xg,yg);
    fi=0;
    PRD=[20 35 45 60];
    for i=PRD
        fi=fi+1;
        periods=unique(Aniso.period);
        ids=find(periods==i);        
        Visos=Aniso.viso(:,ids);
        FPD=Aniso.Fai(:,ids);
        MOA=Aniso.M(:,ids);
        MOA(MOA<0)=0;
        MOA(MOA>5)=5;
        a=MOA.*sin(FPD/180*pi);
        b=MOA.*cos(FPD/180*pi);
        [n,m]=size(xg);
        MOA(Visos<3)=0;
        Visos(Visos<3)=nan;
        for isg=1:n
            for jsg=1:m
                idx=find(x==xg(isg,jsg) & y==yg(isg,jsg));
                if isempty(idx)
                    vGrd(isg,jsg)=nan;
                    ag(isg,jsg)=nan;
                    bg(isg,jsg)=nan;
                else
                    vGrd(isg,jsg)=Visos(idx);
                    ag(isg,jsg)=a(idx);
                    bg(isg,jsg)=b(idx);
                end 
            end
        end
         vGrd=griddata(x,y,Visos,xg,yg);
%         FPDGrd=griddata(x,y,FPD,xg,yg);
%         ag=griddata(x,y,a,xg,yg);
%         bg=griddata(x,y,b,xg,yg);
        
%         subplot(2,2,fi)
        figure
        [~,c]=contourf(lon,lat,VinGrd,200);
        set(c,'edgecolor','none');
        % c.LineColor='w';
        colormap(flipud(jet))
        hold on
        h1=quiver(xg,yg,ag,bg,'k');
        h1.ShowArrowHead = 'off';
        h2=quiver(xg,yg,-ag,-bg,'k');
        h2.ShowArrowHead = 'off';
        hold off
        axis equal
    end