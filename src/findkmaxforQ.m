
function kmax=findkmaxforQ(dat)
maxi = find(dat==max(dat));
diff0=diff(dat); 
ndf=length(diff0);
nmax=0;
pmax = -12345;
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
end
mi=find(ismember(pmax,maxi)==1);
if ~isempty(mi)
    kmax=pmax(mi);
else
    kmax=[];
end