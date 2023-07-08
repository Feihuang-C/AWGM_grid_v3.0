%--------------------------------
%------read dsp for all blocks----------------
%--------------------------------
function dsp=read_dsp(fli)
fid=fopen(fli,'r');
fgets(fid)
fgets(fid)
fgets(fid)
fgets(fid)
fgets(fid)

ibk=0
hh=fscanf(fid,'%g\n',[4])
while(~feof(fid))
    ibk=ibk+1;
    dsp(ibk).ibk=floor(hh(1));
    dsp(ibk).xbk=hh(2);
    dsp(ibk).ybk=hh(3);
    np=floor(hh(4));
    for ir=1:np
        vv=fscanf(fid,'%g\n',[5]);
        pr(ir)=vv(1);
        vg(ir)=vv(2);
        rsl(ir)=vv(3);
        ak(ir)=vv(4);
        bk(ir)=vv(5);
    end
    dsp(ibk).pr=pr;
    dsp(ibk).vg=vg;
    dsp(ibk).rsl=rsl;
    dsp(ibk).ak=ak;
    dsp(ibk).bk=bk;
    hh=fscanf(fid,'%g\n',[4]);
end
fclose(fid);