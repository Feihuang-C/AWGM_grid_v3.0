function writesac0(fli,sacfl,a,b,c,d)
    a(isnan(a)) = -12345.0;
    b(isnan(b)) = -12345;
    for i=1:192
        c0(i)=' ';
    end
    for i=1:24
        c1=cell2mat(c(i));
        nc=length(c1);
        if nc>0
            c0((i-1)*8+1:(i-1)*8+nc)=c1(1:nc);
        end
    end
    c=c0;
    
   
    
    fid = fopen(sacfl,'w+',  'ieee-le');
    counta = fwrite( fid, a, 'float32' );
    countb = fwrite( fid, b, 'int32' );
    countc = fwrite( fid, c, 'char' );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Read SAC file data.
    %
    countd = fwrite( fid, d,'float32' );
    %counta = fwrite( fid, '%f' ,a );
    %countb = fwrite( fid, '%d' ,b );
    %countc = fwrite( fid, '%s' ,c );
    %countd = fwrite( fid, '%f' ,d );
    fclose(fid);
    
   
    