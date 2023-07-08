%%
%% Read dispersion curves %%%
%%
function falt=readfault(fl)
dinv=[];
nd=0;
str=[];
try
fli=cell2mat(fl);
catch
fli=fl;
end
fid=fopen(fli,'r');

while feof(fid)==0
        str=fgetl(fid); 
       
    if ~isempty(strfind(str,'#'));
        continue;
    end
        
    if ~isempty(strfind(str,'>'));
        nd=nd+1;
        ii=0;
%         falt(nd).name=str(2:end);
        falt(nd).name=nd;

    else
        ii=ii+1;
        falt(nd).a(ii,:)=str2num(str);
    end

end
    fclose(fid);
    

