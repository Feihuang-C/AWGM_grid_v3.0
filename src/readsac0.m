% function [A,B,C,D]=readsac0(sacfl)  
%     if strcmp(computer,'GLNX86')
%         fid=fopen(sacfl,'r',  'ieee-le');
%     elseif strcmp(computer,'PCWIN')
%         fid=fopen(sacfl,'r',  'ieee-le');
%     else
%         fid=fopen(sacfl,'r',  'ieee-be');
%     end
%     A = fread( fid, [ 70 1 ], 'float32' );
%     B = fread( fid, [ 40 1 ], 'int32' );
%     C = char(fread ( fid, [ 1 192 ], 'char' ));
%     
%     A(A==-12345.0) = NaN;
%     B(B==-12345) = NaN;
%     C = cellstr(reshape(C,8,24)');
%     C(strmatch('-12345',C)) = {''};
%     if nargout==4
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Read SAC file data.
%         %
%         D = fread( fid, 'float32' );
%     end
%     fclose(fid);





function [A,B,C,D]=readsac0(sacfl)  
    fid=fopen(sacfl,'r',  'ieee-le');
    A = fread( fid, [ 70 1 ], 'float32' );

    dt=A(1);
    if dt<0.000001 | dt>1000
        fclose(fid);
        fid=fopen(sacfl,'r',  'ieee-be');
        A = fread( fid, [ 70 1 ], 'float32' );
    end
    
    B = fread( fid, [ 40 1 ], 'int32' );

    C = char(fread ( fid, [ 1 192 ], 'char' ));
    
    A(A==-12345.0) = NaN;
    B(B==-12345) = NaN;
    C = cellstr(reshape(C,8,24)');
    C(strmatch('-12345',C)) = {''};
    if nargout==4
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Read SAC file data.
        %
        D = fread( fid, 'float32' );
    end
    fclose(fid);