classdef Entry < handle
% Entry Object representing a single entry in a zip archive file.
%
    
    properties(Dependent = true)
        
        MadeByZipVersion    (1,1) char;
        MadeByPlatform      (1,1) char;
        ExtractMinVersion   (1,1) char;
        ExtractMinFeature   (1,1) char;
        CanBeExtracted      (1,1) logical;
        IsEncrypted         (1,1) logical;
        IsStronglyEncrypted (1,1) logical;
        UsesUTF8            (1,1) logical;
        CompressionMethod   (1,:) char;
        FileTimeStamp       (1,:) char;
        CompressedSize      (1,1) uint64
        UncompressedSize    (1,1) uint64
        IsDirectory         (1,1) logucal;
        IsTextFile          (1,1) logical;
        Path                (1,:) char;
        FileName            (1,:) char;

    end
    
    properties(Access = protected)

        % VersionMadeBy            B16(1);
        % VersionNeededToExtract   B16(2);
        % CompressionMethod        B16(3); 
        % LastModFileTime          B16(4);
        % LastModFileDate          B16(5);
        % FileNameLength           B16(6);
        % ExtraFieldLength         B16(7);
        % FileCommentLength        B16(8);
        % DiskNumberStart          B16(9);
        % InternalFileAttributes   B16(10);

        % CRC32                    B32(1);
        % CompressedSize           B32(2);
        % UncompressedSize         B32(3);
        % ExternalFileAttributes   B32(4);
        % Offset                   B32(5);

        % CompessedSize64          B64(1)
        % UncompessedSize64        B64(2)
        % Offset64                 B64(3)
        % FileOffset64             B64(4)
        % DiskNumberStart64        B64(4)

        GeneralPurposeBit_      (1,1) io.GeneralPurposeBit
        RawFileName_            (:,1) uint8
        FileComment_            (:,1) uint8;
        ExtraFields_            (:,1) io.Fields.Field;
        
        B16_                    (10,1) int16;
        B32_                    (5,1) int32;
        B64_                    (5,1) double;
        
        IsDirectory_            (1,1) logical = false;
        FileName_               (1,:) char;
        Path_                   (1,:) char;
    end
    
    properties(Access = private, Constant = true )

    end
    
    methods
        function this = Entry( fid, readLocalHeader )
            
            if nargin
                if nargin == 1
                    readLocalHeader = false;
                end
                
                this.loadFromFile(fid,readLocalHeader);
                    
            end

        end

    end
    
    methods(Access = protected)
        
        function loadFromFile( this, fid, readLocalHeader )
           
            try
                [B16,B32,GPB,rawFileName,extraData,commentField] =...
                    io.internal.parseCentralDirectoryRecord( fid );
            catch ME
                throwAsCaller(ME)
            end
            
            fileName = char(rawFileName');
            
            this.IsDirectory_ = io.Util.isDirectory(fileName);
            
            [this.FileName_,this.Path_] = io.Util.getFileNameFromURL(fileName);
            
            this.RawFileName_ = rawFileName;

            this.FileComment_ = commentField;            

            if ~this.IsDirectory_
                
                neededFlags = [B32(2)==0xFFFFFFFF,B32(3)==0xFFFFFFFF,B32(5)==0xFFFFFFFF,B16(9)<0];

                %fields = org.apache.commons.compress.archivers.zip.ExtraFieldUtils.parse(extraFields',false);

                if any(neededFlags)
                    if isempty(extraData)
                         ME = MException('Entry:BadZip64Extra',...
                            'The Zip64 extended information extra field is required, but not present'); 
                        throwAsCaller(ME);                   
                    end

                    this.ExtraFields_ = io.internal.parseExtraData(extraData);
                    z64ExtraFieldMask = ismember([this.ExtraFields_(:).Id],1);

                    if sum(z64ExtraFieldMask) == 0
                        ME = MException('Entry:BadZip64Extra',...
                            'The Zip64 extended information extra field is required, but not present'); 
                        throwAsCaller(ME);                
                    end

                    z64ExtraField = this.ExtraFields_(z64ExtraFieldMask);

                    [msg,compressedSize,uncompressedSize,offset,diskNumberStart] =...
                        z64ExtraField.update(B32(2),B32(3),B32(5),B16(9));

                    if ~isempty(msg)
                        ME = MException('Entry:BadZip64Extra',msg); 
                        throwAsCaller(ME); 
                    end

                    this.B64_ = [compressedSize;uncompressedSize;offset;offset;diskNumberStart];
                else
                    this.B64_ = double([B32([2,3]);B32(5);B32(5);B16(9)]);
                end

                this.GeneralPurposeBit_ = io.GeneralPurposeBit(GPB);

                this.B16_ = B16;
                this.B32_ = B32;

                if readLocalHeader

                    endOfCD = ftell(fid);

                    fseek(fid,this.Offset_,'bof');

                    status = io.Util.isa(fread(fid,[1,4],'*uint8'),'lfh');

                    if status
                        % Skip over LFH to file name length field
                        fseek(fid,22,'cof');
                        % Read both lengths
                        fileNameLen = fread(fid,1,'uint16');
                        extrasFieldLen = fread(fid,1,'uint16'); % file name length
                        % Read file name
                        fileName = fread(fid,fileNameLen,'*uint8');
                        if fileName ~= this.RawFileName_ 
                            ME = MException('Entry:CorruptFile',...
                                'Data in the central file directory and corresponding local file header does not match\nLikely the zip file is corrupt'); 
                            throwAsCaller(ME);  
                        end

                        % And extra fields
                        this.ExtraFields_ = [this.ExtraFields_,...
                            io.internal.parseExtraData(fread(fid,extrasFieldLen,'*uint8'))];

                        this.FileOffset_ = uint64(ftell(fid));

                    end

                    this.B64_(4) = uint64(ftell(fid));

                    fseek(fid,endOfCD,'bof');

                end
            
            end
   
        end
    end
   
    methods
       
        function fullFileName = getFullFileName( this )
            if this.GeneralPurposeBit_.useUTF8ForNames
                fullFileName = native2unicode(char(this.RawFileName_'),'UTF-8');
            else
                fullFileName = char(this.RawFileName_');
            end
        end
        
        function fileName = getFileName( this )
            if this.GeneralPurposeBit_.useUTF8ForNames
                fileName = native2unicode(this.FileName_,'UTF-8');
            else
                fileName = this.FileName_;
            end         
        end
        
        function path = getPath( this )
            if this.GeneralPurposeBit_.useUTF8ForNames
                path = native2unicode(this.Path_,'UTF-8');
            else
                path = this.Path_;
            end        
        end   
        
        function offset = getFileOffset( this )
            offset = this.FileOffset_;
        end

        function gpb = getGPB( this )
           gpb = this.GeneralPurposeBit_; 
        end
        
        function code = getMethod( this )            
            code = uint16(this.B16_(3));
        end
        
        function code = getPlatformCode( this )
            bytes = typecast(this.B16_(1),'uint8');
            code = bytes(2);
        end
        
        function mode = getFileMode( this )
            
            switch this.getPlatformCode
                case 0
                    mode = this.ExternalFileAttributes_;
                case 3 
                    mode = bitand(bitshift(this.ExternalFileAttributes_,-16,'uint32'),l);
                otherwise
                    mode = 0;
            end
            
        end
        
        function data = getRawExtraData( this )
            data = this.ExtraFields_;
        end
        
        function comments = getFileComments( this )
            if isempty(this.FileComment_)
                comments = '';
            else
                comments = char(this.FileComment_');
            end
        end
        
    end
    
    methods
       
        function value = get.MadeByZipVersion( this )
            bytes = typecast(this.B16_(1),'uint8');
            value = io.Util.getVersion(bytes(1));
        end
        
        function value = get.MadeByPlatform( this )
            bytes = typecast(this.B16_(1),'uint8');
            value = io.Util.getPlatform(bytes(2));
        end   
        
        function value = get.ExtractMinVersion( this )
            bytes = typecast(this.B16_(2),'uint8');
            value = io.Util.getVersion(bytes(1));
        end   
        
        function value = get.ExtractMinFeature( this )
            bytes = typecast(this.B16_(2),'uint8');        
            value = io.Util.getMinimumFeature(bytes(2));
        end      
            
        function value = get.IsEncrypted( this )
            value = this.GeneralPurposeBit_.useEncyrption;
        end
        
        function value = get.IsStronglyEncrypted( this )
            value = this.GeneralPurposeBit_.useStrongEncyrption;
        end        

        function value = get.UsesUTF8( this )
            value = this.GeneralPurposeBit_.useUTF8ForNames;
        end 
        
        function value = get.CompressionMethod( this )
            value = io.Util.getCompressionMethod(this.B16_(3));
        end 

        function value = get.CompressedSize( this )
            value = this.B64_(1);
        end

        function value = get.UncompressedSize( this )
            value = this.B64_(2);
        end
        
        function fileName = get.FileName( this )
            fileName = this.getFileName();            
        end  
        
        function path = get.Path( this )
            path = this.getPath();            
        end 
        
        function value = get.IsDirectory( this )
            value = this.IsDirectory_;
        end         
        
        function boolean = get.IsTextFile( this )
           boolean = logical(bitget(this.B16_(10),1,'int16')); 
        end
        
        function boolean = get.CanBeExtracted( this )
            supportedModes = io.Util.getSupportedModes();
            isExtractableByMode = ismember(this.CompressionMethod_,supportedModes);
            
            boolean = isExtractableByMode && this.GeneralPurposeBit_.isExtractable();            
        end        

    end
end



