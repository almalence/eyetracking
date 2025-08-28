function [B16,B32,GPB,rawFileName,extraFields] =...
        parseLocalFileHeader( fid )

    B16 = zeros([6,1],'int16');
    B32 = zeros([3,1],'int32');

    B16(1) = fread(fid,1,'*int16');
    GPB = io.GeneralPurposeBit(fread(fid,16,'ubit1=>uint8'));
    B16(2:4) = fread(fid,3,'*int16');
    B32(1:3) = fread(fid,3,'*int32');
    B16(5:6) = fread(fid,2,'*int16');            

    rawFileName = fread(fid,[1,B16(5)],'*uint8');

    if B16(6) > 0
        extraFields = io.internal.parseExtraData(fread(fid,B16(6),'*uint8'));
    else
        extraFields = [];
    end


end  

