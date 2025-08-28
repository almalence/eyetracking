function [noEntries, extraData, comments] = positionAtCentralDirectoryRecord( fid )
% positionAtCentralDirectoryRecord - Attempts to position the file
% pointer at the first central directory record by locating and then 
% parsing the end of central directory record or Zip64 end of 
% central directory locator.
%
% inspired by:  
% org.apache.commons.compress.archivers.zip.ZipFile
%
% Reference:
% https://pkware.cachefly.net/webdocs/casestudies/APPNOTE.TXT

    extraData = [];
    comments = [];

    % Set pointer to eof to get file length
    fseek(fid,0,'eof');    
    fileLength = ftell(fid);
    % Zip files cannot be smaller than the length of the end of centeral
    % directory record (22 btes )
    if fileLength < io.Util.MIN_LENGTH
        ME = MException('getZipContents:NotAZip','Archive is not a ZIP archive, file length is less than 22 bytes');
        throwAsCaller(ME);
    end
    % Start search at EOCD length(22 bytes) from eof
    eocdPos = -double(io.Util.EOCD_CHUNK_SIZE);
    fseek(fid,eocdPos,'eof');
    % Check for EOCD signature
    foundEOCDSig = io.Util.isa(fread(fid,4,'*uint8'),'eocd');
    % Search while EOCD signature is not found and we have not searched
    % farther than the maximum (supposedly) length of the EOCD
    % block
    % ZIP_MAGIC_NUMBER = 65535 bytes.
    % See reference 4.4.12
    % TODO: Maybe faster just to read in MAGIC NUMBER of bytes and
    % search for EOCD signature in that byte array, would depend on
    % length of comment field
    maxSearchLength = -min(fileLength,io.Util.MAX_LENGTH);
    while (~foundEOCDSig) && (eocdPos > maxSearchLength)
        eocdPos = eocdPos - 1;
        fseek(fid,eocdPos,'eof');
        foundEOCDSig = io.Util.isa( fread(fid,4,'*uint8'), 'EOCD');               
    end
    % Error if EOCD still not found
    if ~foundEOCDSig        
        ME = MException('getZipContents:NotAZip','Archive is not a ZIP archive or is corrupt.\nEnd of centeral directory signature could not be found');
        throwAsCaller(ME);
    end
    % Reset pointer back to start of EOCD and then read EOCD as a
    % block
    fseek(fid,eocdPos,'eof');
    [~,bytes] = io.Util.readChunk( fid, 'eocd'); 
    offset = typecast(bytes(17:20),'uint32');
    noEntries = typecast(bytes(9:10),'uint16'); 
    commentLen = typecast(bytes(21:22),'uint16');
    if commentLen > 0
        comments = fread(fid,[1,commentLen],'*char');
    end
    % Move pointer towards bof to search for ZIP64EOCD locator
    % signature.
    % This assumes that the EOCD is preceded by the Zip64EOCDL
    % chunk as does apache.commons.compress.zip.ZipFile and is
    % hinted at in the reference
    fseek(fid,eocdPos-double(io.Util.Z64EOCDL_CHUNK_SIZE),'eof');
    % Get ZIP64EOCD locator chunk
    [foundZ64EOCDLSig,Z64EOCDLBytes] = io.Util.readChunk(fid,'z64eocdl');
    % If  Zip64 EOCDL is found (ver 1/2 (?) format)
    if foundZ64EOCDLSig
        % Read Zip64 end of central directory locator and set 
        % pointer to Zip64 end central directory record, then parse
        % the record
        %
        offsetToZ64EOCDR = typecast(Z64EOCDLBytes(9:16),'uint64');
        % Set pointer to start of Zip64 end of centeral directory record
        fseek(fid,offsetToZ64EOCDR,'bof'); 
        % Parse Zip64 end of central directory locator, find offset
        % to central directory record
        [offset,noEntries,extraData] = io.internal.readZip64EndOfCentralDirectoryRecord( fid );
    end
    % Set pointer to start of central directory
    fseek(fid,offset,'bof'); 
end  
