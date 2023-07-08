
%-------------------new num2str-----------------------
function str=num2str2(val,nb,ne);
s0=num2str(val);
if round(val)==val
    se(1:ne)='0';
    if(ne==0)       
    else
        s0=[s0 '.' se];
    end
else
    ii=strfind(s0,'.');
    nn=length(s0);
    ne0=nn-ii;
    if ne0<ne
        s0(nn+1:nn+ne-ne0)='0';
    elseif ne0>ne
        s0(nn-ne0+ne+1:nn)='';
    end
end
if length(s0)<nb
  
    str(1:nb-length(s0))='0';
    str(nb-length(s0)+1:nb)=s0;
elseif length(s0)>nb
   
    str(1:nb)=s0(1:nb);
else
    str=s0;
end

%-------------------test only-----------------------
function test
dd=2.011;
str=num2str2(dd,5,2);

