function [offsetToCD,noEntriesOnThisDisk,extensibleData] =...
        readZip64EndOfCentralDirectoryRecord( fid )
% readZip64EndOfCentralDirectoryRecord - Parse Zip64 end of central
% directory record for 8 byte offset to centeral directory and 
% number of entries, return extnesible data if present

    [status,bytes ] = io.Util.readChunk( fid, 'z64eocdr' );

    if status

        sizeOf = typecast(bytes(5:12),'uint64');
        % verMadeBy = typecast(bytes(12:14),'uint16');
        % verToExtract = typecast(bytes(15:16),'uint16');
        % noOfDisks = typecast(bytes(17:20),'uint32');
        % noDiskWSCD = typecast(bytes(21:24),'uint32');
        noEntriesOnThisDisk = typecast(bytes(25:32),'uint64');
        % noEntriesTotal = typecast(bytes(33:40),'uint64');
        % szCD = typecast(bytes(41:48),'uint64');
        offsetToCD = typecast(bytes(49:56),'uint64');

        sizeOfExtensibleData = sizeOf - (uint64(io.Util.Z64EOCDR_CHUNK_SIZE) + uint64(12));

        extensibleData = [];

        if sizeOfExtensibleData > 0
            extensibleData = fread(fid,sizeOfExtensibleData,'*uint8');
        end


    else
        % TODO: See if we can recover from this error
        ME = MException('readZip64EndOfCentralDirectoryRecord:BadSignature',...
            'Could not find the Zip64 end of central directory record'); 
        throwAsCaller(ME);                 
    end         

end  

