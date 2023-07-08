
function vmod=smoothmod(vmodin,drmax)
ndp=length(vmodin(1).dp);
lo=[vmodin.xbk];
la=[vmodin.ybk];
visos=cell2mat({vmodin.viso}');
as=cell2mat({vmodin.a}');
bs=cell2mat({vmodin.b}');
FPDs=cell2mat({vmodin.FPD}');
MOAs=cell2mat({vmodin.MOA}');
vmod=vmodin;
for i=1:ndp
    disp(['dp:  ' num2str(i)])
        viso=visos(:,i);
        a=as(:,i);
        b=bs(:,i);
        for k=1:length(a)
          loi=lo(k);
          lai=la(k);
          dist=distance(lai,loi,la,lo);
          id=find(dist<=drmax);
          vi=mean(viso(id));
          ai=mean(a(id));
          bi=mean(b(id));
          azm=0:179;
          vazms=fun1([vi,ai,bi],azm/180*pi);
          FPDi=azm(vazms==max(vazms));
          FPDi=FPDi(1);
          MOAi=(max(vazms)-min(vazms))/vi*100;          
          vmod(k).viso(i)=vi;
          vmod(k).a(i)=ai;
          vmod(k).b(i)=bi;
          vmod(k).FPD(i)=FPDi;
          vmod(k).MOA(i)=MOAi;
          vmod(k).vs(i,:)=fun1([vi,ai,bi],(0:10:170)/180*pi);          
          vmod(k).azm(i,:)=0:10:170;          
        end
        
    
end
  
function fy1=fun1(abc,fx)
fy1=abc(1)+abc(2)*cos(2*fx)+abc(3)*sin(2*fx);