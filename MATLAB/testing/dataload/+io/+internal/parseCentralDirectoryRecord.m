function [B16,B32,GPB,rawFileName,extraFields,commentField] =...
        parseCentralDirectoryRecord( fid )

    B16 = zeros([10,1]      ,'int16');
    B32 = zeros([5,1]       ,'uint32');

    B16(1:2) = fread(fid,2  ,'*int16');
    GPB = fread(fid,16      ,'ubit1=>uint8');
    B16(3:5) = fread(fid,3  ,'*int16');
    B32(1:3) = fread(fid,3  ,'*uint32');
    B16(6:10) = fread(fid,5 ,'*int16');
    B32(4:5) = fread(fid,2  ,'*uint32');               

    rawFileName = fread(fid,B16(6),'*uint8');

    if B16(7) > 0
        extraFields = fread(fid,B16(7),'*uint8');
    else
        extraFields = [];
    end

    if B16(8) ~= 0
        commentField = fread(fid,B16(8),'*char');
    else
        commentField = [];
    end 

end  

