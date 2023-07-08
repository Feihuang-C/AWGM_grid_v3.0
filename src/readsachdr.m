function [A,B,C]=readsachead(sacfl)  
    if strcmp(computer,'GLNX86')
        fid=fopen(sacfl,'r',  'ieee-le');
    else
        fid=fopen(sacfl,'r',  'ieee-be');
    end
    A = fread( fid, [ 70 1 ], 'float32' );
    B = fread( fid, [ 40 1 ], 'int32' );
    C = char(fread ( fid, [ 1 192 ], 'char' ));
    
    fclose(fid);
