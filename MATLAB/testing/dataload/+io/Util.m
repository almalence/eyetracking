classdef Util

    properties(Constant = true)
        
        
        ZIP64_EOCD_LOC      = uint8([80;75;6;7]);
        ZIP64_EOCD_REC      = uint8([80;75;6;6]);
        CFH                 = uint8([80;75;1;2]);
        LFH                 = uint8([80;75;3;4]);
        EOCD                = uint8([80;75;5;6]);        
        DATA_DESCRIPTOR     = uint8([80;75;7;8]);
        DIGITAL             = uint8([80;75;5;5]);
        AED                 = uint8([80;75;6;8]);

        Z64_EXT_SIZE        = [8;8;8;4]; 
        Z64_EXT_FIELD       = uint8([1;0]);
        
        SHORT               = uint8(2);
        WORD                = uint8(4);
        DWORD               = uint8(8);
        MAGIC_NUMBER        = uint32(65535);
        ZIP64_MAGIC_SHORT   = uint32(65535);        
        ZIP64_EOCDL_LENGTH  = uint32(20);
        MIN_LENGTH          = uint32(22);
        MAX_LENGTH          = 65535; 
        
        EOCD_CHUNK_SIZE     = uint32(22);
        CFH_CHUNK_SIZE      = uint32(46);
        Z64EOCDR_CHUNK_SIZE = uint32(56);
        Z64EOCDL_CHUNK_SIZE = uint32(20);
        LFH_CHUNK_SIZE      = uint32(30);
        DIGITAL_CHUNK_SIZE  = uint32(6);
    end
    
    methods(Static = true)
%--------------------------------------------------------------------------         
        function bytes = hex2LongBytes( hexStr )
        % hex2bytes returns a row vector of uint8 bytes from a hex string
%             bytes = uint8(org.apache.commons.compress.archivers.zip.ZipLong.getBytes(hex2dec(hexStr))');
            bytes = typecast(uint32(hex2dec(hexStr)),'uint8');
        end  
%--------------------------------------------------------------------------       
        function bytes = hex2ShortBytes( hexStr )
        % hex2bytes returns a row vector of uint8 bytes from a hex string
%             bytes = uint8(org.apache.commons.compress.archivers.zip.ZipShort.getBytes(hex2dec(hexStr))');
            bytes = typecast(uint16(hex2dec(hexStr)),'uint8');
        end        
%--------------------------------------------------------------------------        
        function [status, position] = isa( bytes, strId )

            position = [];
            noBytes = numel(bytes);
            if noBytes == 4
                status = isequal(bytes,io.Util.get(strId));
            elseif noBytes > 4
                position = strfind(bytes,io.Util.get(strId));
                status = ~isempty(position);
                % First occurence
%                 position = position(1);
            else
                status = false;
            end
        end
%--------------------------------------------------------------------------        
        function sig = get( strId )

            switch lower(strId)

                case 'eocd' % End of central directory signature
                    sig = io.Util.EOCD;
                case 'z64eocdl' % zip64 end of central directory locator signature
                    sig = io.Util.ZIP64_EOCD_LOC;
                case 'z64eocdr' % zip64 end of centeral directory record signature
                    sig = io.Util.ZIP64_EOCD_REC;
                case 'digital'
                    sig = io.Util.DIGITAL;
                case 'cfh' % Central file header signature
                    sig = io.Util.CFH;
                case 'aedr' % Archive extra data signature  
                    sig = io.Util.AED;
                case 'lfh' % Local file header signature 
                    sig = io.Util.LFH;
                case 'z64ext'
                    sig = io.Util.Z64_EXT_FIELD; 
                otherwise
                    sig = uint8([255,255,255,255]);
            end
        end 
%--------------------------------------------------------------------------        
        function [ status, bytes ] = readChunk( fid, chunkType )

            switch lower(chunkType)
                case 'eocd'
                    chunkSize = io.Util.EOCD_CHUNK_SIZE;
                case 'z64eocdr'                    
                    chunkSize = io.Util.Z64EOCDR_CHUNK_SIZE;
                case 'z64eocdl'
                    chunkSize = io.Util.Z64EOCDL_CHUNK_SIZE;
                case 'cfh'                    
                    chunkSize = io.Util.CFH_CHUNK_SIZE;                    
                case 'lfh'
                    chunkSize = io.Util.LFH_CHUNK_SIZE;
                case 'digital'
                    chunkSize = io.Util.DIGITAL_CHUNK_SIZE;
                otherwise
                    error('a:b','c');
            end
            
            try
                bytes = fread(fid,chunkSize,'*uint8');
                status = io.Util.isa(bytes(1:4),chunkType);
            catch
                bytes = [];
                status = false;
            end

        end
%--------------------------------------------------------------------------        
        function URL = constructURL( fullFileName, archiveFileName )

            if startsWith(fullFileName,'/')
                fullFileName = fullFileName(2:end);
            end
            archiveFileName = strrep(archiveFileName,'\','/');
            URL = ['jar:file:/',archiveFileName,'!/',fullFileName];
        end
%--------------------------------------------------------------------------        
        function [fileName, path] = getFileNameFromURL( URL )
            
            if ~contains(URL,'/')
                fileName = URL;
                path = '';
                return
            end
            
            pos = find(URL == '/',1,'last') + 1;
            fileName = URL(pos:end);
            path = URL(1:pos-1);
        end   
%--------------------------------------------------------------------------        
        function tf = isDirectory( fileName )
            tf = endsWith(fileName,'/');
        end
%--------------------------------------------------------------------------        
        function n = validateNumber( n, fcnName )

            try
                validateattributes(n,{'numeric'},...
                    {'scalar','positive','integer','finite'},fcnName,'n',1)
            catch ME
                thowAsCaller(ME)
            end
        end  
%--------------------------------------------------------------------------        
        function handleIOExceptions( ME )
            % TODOD: Meaningful error messages
            if isa(ME,'java.io.IOException')
                cause = char(ME.getCause());
                if isempty(cause)
                    cause = 'Unkown cause';
                end
                message = char(ME.getMessage);
                if isempty(message)
                    message = 'No details';
                end
                
                ME = MException('IOException','%s\n%s',cause,message);
            end
            throwAsCaller(ME);
        end
%--------------------------------------------------------------------------        
        function time = dosTimeToJavaTime( bytes )
            import org.apache.commons.compress.archivers.zip.*

            value = ZipLong.getValue(bytes);
            time = ZipUtil.dosToJavaTime(value);

        end  
%--------------------------------------------------------------------------        
        function version = getVersion( code )

            val = double(code);

            minver = mod(val,10);
            
            maxver = val/10;

            if minver == 0
                if maxver == 0
                   version = 'Not defined'; 
                else
                    version = sprintf('%2.1f',val/10);
                end
            else
                version = sprintf('%2.1f.%2.1f',val/10,minver);
            end
        end
%--------------------------------------------------------------------------
        function platform = getPlatform( code )

            switch code
                case 0  ;platform = 'FAT/VFAT/FAT32';
                case 1  ;platform = 'Amiga';
                case 2  ;platform = 'OpenVMS';                
                case 3  ;platform = 'UNIX';
                case 4  ;platform = 'VM/CMS';                 
                case 5  ;platform = 'Atari ST';
                case 6  ;platform = 'OS/2 H.P.F.S.';                 
                case 7  ;platform = 'Macintosh';
                case 8  ;platform = 'Z-System';                 
                case 9  ;platform = 'CP/M';
                case 10 ;platform = 'Windows NTFS';                 
                case 11 ;platform = 'MVS (OS/390 - Z/OS)';
                case 12 ;platform = 'VSE';                 
                case 13 ;platform = 'Acorn Risc';
                case 14 ;platform = 'VFAT';                 
                case 15 ;platform = 'alternate MVS';
                case 16 ;platform = 'BeOS';                 
                case 17 ;platform = 'Tandem';
                case 18 ;platform = 'OS/400';                 
                case 19 ;platform = 'OS X (Darwin)';
                otherwise ;platform = 'undefined';
            end

        end
%--------------------------------------------------------------------------
        function feature = getMinimumFeature( code )

            switch code
                case 0 ; feature = 'None required';
                case 1 ; feature = 'Default';
                case 11; feature = 'Volume label';
                case 20; feature = 'Deflate compression';
                case 21; feature = 'Deflate64';
                case 25; feature = 'PKWARE DCL Implode';
                case 27; feature = 'Patch data set';
                case 45; feature = 'ZIP64 format extension';
                case 46; feature = 'BZIP2 Compression';
                case 50; feature = 'DES, 3DES, original RC2, or RC4';
                case 51; feature = 'AES or corrected RC2 encryption';
                case 52; feature = 'Corrected RC2-64 encryption';
                case 61; feature = 'non-OAEP key wrapping';
                case 62; feature = 'Central directory encryption';
                case 63; feature = 'LZMA, PPMd+, Blowfish or Twofish';
                otherwise ;feature = 'undefined';
            end

        end
%--------------------------------------------------------------------------
        function method = getCompressionMethod( code )

            switch code
                case 0; method = 'No compression';
                case 1; method = 'Shrunk';
                case 2; method = 'Reduced with compression factor 1';
                case 3; method = 'Reduced with compression factor 2';
                case 4; method = 'Reduced with compression factor 3'; 
                case 5; method = 'Reduced with compression factor 4';
                case 6; method = 'Imploded';
                case 7; method = 'Reserved for Tokenizing compression algorithm';
                case 8; method = 'Deflated';
                case 9; method = 'Defalte64';
                case 10; method = 'PKWARE Data Compression Library Imploding (old IBM TERSE)';
                case 12; method = 'BZIP2';
                case 14; method = 'LZMA';
                case 16; method = 'IBM z/OS CMPSC Compression';
                case 18; method = 'File is compressed using IBM TERSE';
                case 19; method = 'IBM LZ77 z Architecture';
                case 20; method = 'Zstandard (zstd) Compression';
                case 96; method = 'JPEG variant';
                case 97; method = 'WavPack compressed data';
                case 98; method = 'PPMd version I, Rev 1';
                case 99; method = 'AE-x encryption marker';
                otherwise; method = 'Undefined';
            end

        end        
%--------------------------------------------------------------------------     
        function exts = getImageIOFileExtensions()

            % Create a cell array of image extension recognised by ImagIO
            rdrSuffixes = javax.imageio.ImageIO.getReaderFileSuffixes();
            noSuffixes = numel(rdrSuffixes);
            exts = cell(1,noSuffixes);
            for itr = 1:noSuffixes
                exts{itr} = ['.',char(rdrSuffixes(itr))]; 
            end
        end    
%--------------------------------------------------------------------------   
        function ver = getCommonCompressVersion()

            javaPaths = javaclasspath('-all');
            compressJarFiles = javaPaths(contains(javaPaths,'commons-compress'));
            if isempty(compressJarFiles)
                ver = 'Unknown';
                return
            end
            numJarFiles = numel(compressJarFiles);
            ver = cell(1,numJarFiles);
            k = 1;
            for itr = 1:numJarFiles
                compressJarFile = compressJarFiles{itr};
                if exist(compressJarFile,'file') == 2
                    compressJarFile = strrep(compressJarFile,'\','/');
                    URLstr = ['jar:file:/',compressJarFile,'!/'];
                    URL = java.net.URL(URLstr);
                    jarURLconnection = URL.openConnection;
                    manifest = jarURLconnection.getManifest();
                    attribs = manifest.getMainAttributes();
                    ver{k} = char(attribs.getValue('Implementation-Version'));
                    k = k + 1;
  
                end
            end
            
            if numJarFiles == 1; ver = ver{1}; end
            
            
        end  
%--------------------------------------------------------------------------  
        function modes = getSupportedModes()
            persistent savedModes
            
            if isempty(savedModes)
                ver = io.Util.getCommonCompressVersion();
                switch ver
                    case '1.8.1'
                        savedModes = uint16([0,1,6,8]);
                    case '1.20'
                        savedModes = uint16([0,1,6,8,9,12]);
                    otherwise
                        savedModes = [];
                end
            end
            
            modes = savedModes;
            
        end 
%-------------------------------------------------------------------------- 
        function [ zipFile, fid, cleanUpObj ] = validateZipFile( zipFile )
            
            % Open the file for reading
            [fid,msg] = fopen(zipFile,'r');            

            if fid < 3
               ME = MException('getZipContents:badFile','Cannot read the file %s\nfoopen returned : %s',zipFile,msg);
               throwAsCaller(ME);
            end
            zipFile = fopen(fid);
            cleanUpObj = onCleanup(@()fclose(fid));
            % Check first 4 bytes of file for a 'lfh' signature. 
            if ~io.Util.isa( fread(fid,4,'*uint8'), 'lfh' )
               ME = MException('getZipContents:badLFH',...
                   'The file, %s,\ndoes not start with a local file header signature\nIt may be a self-exctracting archive file or an empty zip file',zipFile); 
               throwAsCaller(ME);
            end
        end
%-------------------------------------------------------------------------- 
        function file = uigetzipfile( path )
            
            if nargin == 0
                path = pwd;
            end
            
            filter = {'*.zip','Zip File (*.zip)';...
                      '*.jar;*.war;*.ear;*.rar','Java Files (*.jar,*.war,*.rar,*.rar)';...
                      '*.docx;*.docm','Office Open (*.docx,*.docm)';...
                      '*.*',  'All Files (*.*)'};
                  
            [filename,path] = uigetfile(filter,'Select zip file',path);
            
            if isequal(filename,0)
                file = [];
                return
            end
            
            file = fullfile(path,filename);
        end
    end    
end

