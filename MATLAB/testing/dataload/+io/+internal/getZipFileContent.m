function [entries, extraData, digitalSignature, zipFile, comments ] =...
        getZipFileContent( zipFile, readLocalFileHeader )

    [ zipFile, fid, cleanpObj ] = io.Util.validateZipFile( zipFile ); %#ok<ASGLU>
    
    [noEntries,extraData,comments] = io.internal.positionAtCentralDirectoryRecord( fid );

    if ~isempty(noEntries)
        [entries,digitalSignature] =...
            io.internal.populateFromCentralDirectory( fid, noEntries, readLocalFileHeader );
    else
        % TODO: Better explaination
       ME = MException('getZipContents:NoEntries',...
           'No entries found'); 
       throwAsCaller(ME);            
    end
end 

