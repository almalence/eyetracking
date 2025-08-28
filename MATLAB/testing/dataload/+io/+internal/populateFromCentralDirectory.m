function [entries,digitalSignature] = populateFromCentralDirectory( fid, noEntries, readLocalFileHeader )
% populateFromCentralDirectory - Parses central directory records
% and retuen an array of io.Entry objects

    %Pre-allocation
    entries(noEntries) = io.Entry();
    digitalSignature = [];

    CFH = io.Util.CFH;

    doLoop = true;

    k = 1;

    while doLoop
        % Verify and returns a byte array of the fixed fields
        doLoop = isequal(fread(fid,4,'*uint8'),CFH);

        if doLoop
        % The constructor of the io.Entry objects parses fixed
        % fields and reads any variable length fields so that the file
        % pointer should be at the beginning of the next record
            entries(k) = io.Entry( fid, readLocalFileHeader );

            k = k + 1;
        else
        % We assume doLoop is false because all the central directory
        % entries have been read and then attempt to read the digital 
        % signature block, if present.  

            [ status, bytes ] = io.Util.readChunk( fid, 'digital' );

            if status
                digitalSigSz = typecast(bytes(5,6),'uint16');
                digitalSignature = fread(fid,digitalSigSz,'*uint8');                        
            end
        end

    end
end 
