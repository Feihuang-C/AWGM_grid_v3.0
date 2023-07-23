%% make dispersion files for 3D inversion
%% Cao Feihaung et al. 2023
function af_step5_WG2surface96

clc
clear
fs = '/';
indir='./aniso/';
infile = [indir fs 'Aniso360azmbin10onephs0_allst_smr.mat']; % input file
outpath = './inv2depth2/';
is_parfor=1;
prdLim = [10,80];    % period range of DC
azms = 0:10:170;  % =0,output isotropy dispersion curves(DC); 
                     % =0:10:170 or =0:20:350, output azimuth dependent DC
if ~exist(outpath,'dir')
    mkdir(outpath)
end

temp = load(infile);
Aniso = temp.Aniso;
NANid = isnan(Aniso.a);
Aniso.Fai(NANid) = 0;
Aniso.M(NANid) = 0;
Aniso.stdM(NANid) = 0;
Aniso.a(NANid) = 0;
Aniso.b(NANid) = 0;
Aniso.stdFai(NANid) = 0;
parmat='./wg2srf96_par.mat';
save(parmat)

if is_parfor==1
    parfor blcki = 1:length(Aniso.st(:,1))
        parWG2sruface96(blcki,parmat)
    end
else
    for blcki = 1:length(Aniso.st(:,1))
        parWG2sruface96(blcki,parmat)
    end
end
copyfile('../src/job_submit.sh',outpath)
copyfile('../src/mpi_job_submit.sh',outpath)
copyfile('../src/mod.d2',outpath)
copyfile('../src/rfnt',outpath)
copyfile('../src/stacksac2.0AJWA',outpath)
copyfile('../src/DispInv.sh',outpath)
copyfile('../src/run_shell_mpi.f90',outpath)
copyfile('../src/jointinversion.sh',outpath)
copyfile('../src/Extract_3Daniso/',[outpath '/Extract_3Daniso/'])
delete(parmat)


function parWG2sruface96(blcki,parmat)
disp(blcki)
load(parmat) 
blockID = num2str2(blcki,5,0);
blckF = [outpath blockID fs];
period=Aniso.period;
Viso=Aniso.viso(blcki,:);
dv=Aniso.stdv(blcki,:);
a=Aniso.a(blcki,:);
b=Aniso.b(blcki,:);
delID=isnan(Viso)|Viso==0|period>prdLim(2)|period<prdLim(1);
period(delID)=[];
Viso(delID)=[];
dv(delID)=[];
a(delID)=[];
b(delID)=[];
if isempty(Viso)
    return
end

if ~exist(blckF,'dir')
    mkdir(blckF)
end

fidlc = fopen([blckF 'xyloc.txt'],'w');
fprintf(fidlc,'%4.0f %7.3f %7.3f',str2double(blockID),Aniso.st(blcki,2),Aniso.st(blcki,1));
fclose(fidlc);

if length(azms)==1
    fid = fopen([blckF 'bk.' blockID '.deg.iso.dsp'],'w');
    for pdi=1:length(period)
        fprintf(fid,'SURF96 R C X 0 %7.2f %7.3f %7.3f %7.3f \n',period(pdi),Viso(pdi),dv(pdi),6);
    end
    fclose(fid);
else
    for azi = 1:length(azms)
        V = Viso+a.*cos(2*(azms(azi)/180*pi))+b.*sin(2*(azms(azi)/180*pi));
        fid = fopen([blckF 'bk.' blockID '.deg.' num2str2(azms(azi),3,0),'.dsp'],'w');
        for pdi=1:length(period)
%             fprintf(fid,'SURF96 R C X 0 %7.2f %7.3f %7.3f %7.3f \n',period(pdi),V(pdi),dv(pdi),6);
            fprintf(fid,'SURF96 R C X 0 %7.2f %7.3f %7.3f %7.3f \n',period(pdi),V(pdi),0,6);
        end
        fclose(fid);
    end
end









