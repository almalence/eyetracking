classdef DataReader < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Dependent = true)
        
        FileName (1,:) char
        
        FileLength (1,1) double
       
        BytesRead (1,1) double
    end
    
    properties(Access = protected)
       
        inputStream
        
        Entry_ (1,1) io.Entry
        
        BytesRead_ (1,1) double = 0;
        
        BytesRemaining_ (1,1) double
        
        FileLength_ (1,1) double = 0;
        
        FileName_ (1,:) char = '';
    end
    
    
    methods
        function this = DataReader( input, entry, ~ )
            
            narginchk(0,3);
            
            if isequal(nargin,0)
                
                [input,this.FileLength_,this.FileName_] = loadFromFile();
                
                if isempty(input)
                    return
                end
                       
            elseif isjava(input)
                
                if ~isa(input,'java.io.InputStream')

                end
            elseif ischar(input)
                
                [fid,errmsg] = fopen(input);

                if fid < 3
                    error('listArchiveEntries:badFile',...
                    'Unable to open %s , fopen returns %s',zipFileName,errmsg)        
                end

                fileName = fopen(fid);

                fclose(fid); 
                
                file = java.io.File(fileName);
                
                this.FileLength_ = file.length();
                
                input = org.apache.commons.io.FileUtils.openInputStream(file);
                
            elseif isempty(input) 
                
                [input,this.FileLength_,this.FileName_] = loadFromFile();
                
                if isempty(input)
                    return
                end
            end
            
            this.inputStream = input;

            if (nargin > 1)
                this.Entry_ = entry;
            end
            
            this.BytesRemaining_ = this.FileLength;

        end
        
        function val = get.BytesRead( this )
            val = this.FileLength - this.BytesRemaining_;
        end
        
        function val = get.FileLength( this )
            
            if isequal(this.Entry_.UncompressedSize,0)
                val = this.FileLength_;
            else
                val = this.Entry_.UncompressedSize;
            end
        end 
        
        function val = get.FileName( this )
            
            if isempty(this.Entry_.FileName)
                val = this.FileName_;
            else
                val = this.Entry_.FileName;
            end
        end        

    end
    
    methods(Access = public)
        
        function val = readBoolean( this, n )
            
            val = logical(this.convert( n, 1, 'uint8' ));
        end
        
        function val = readByte( this, n )
            
            val = this.convert( n, 1, 'uint8' );
        end 
        
        function val = readSignedByte( this, n )
            
            val = this.convert( n, 1, 'int8' );
        end        
        
        function val = readChar( this, n )
            
            byte = this.convert( n, 1, 'uint8' );
            
            val = char(byte);
            
            if isvector(val)
                val = val(:)';
            end
        end
        
        function val = readDouble( this, n )

            val = this.convert( n, 8, 'double' );
        end
        
        function val = readSingle( this, n )
            
            val = this.convert( n, 4, 'single' );
        end
        
        function val = readUint16( this, n )
            
            val = this.convert( n, 2, 'uint16' );
        end
        
        function val = readInt16( this, n )
            
            val = this.convert( n, 2, 'int16' );
        end        
        
        function val = readUint32( this, n )
            
            val = this.convert( n, 4, 'uint32');
        end        
        
        function val = readInt32( this, n )
            
            val = this.convert( n, 4, 'int32' );
        end        
        
        function val = readUint64( this, n )
            val = this.convert( n, 8, 'uint64' );
        end
        
        function val = readInt64( this, n )
            
            val = this.convert( n, 8, 'int64' );
        end        
        
        function skipped = skipBytes( this, n )
            
            skipped = [];
            if nargin == 1
                n = 1;
            end
            validateattributes(n,{'numeric'},{'scalar','positive','finite'},...
                'skipBytes');
            try
                if this.isAvailable(n)
                    
                    skipped = org.apache.commons.io.IOUtils.skip(this.inputStream,n);                    
                else
                    
                    this.throwOverRequesrError(n);                    
                end
            catch
                
                this.throwOverRequesrError(n);                
            end
                    
          
        end
        
        function tf = isAvailable( this, bytes )
            
            tf = this.BytesRemaining_ >= bytes;
        end
        
        
    end
    
    methods(Access = protected)
        
        function out = convert( this, n, multi, type )

            if nargin == 1
                n = 1;
                sz = [];
            else
                [n,sz,ME] = this.validateDimension( n );
                if ~isempty(ME)
                    throwAsCaller(ME);
                end
            end
            
            n = n * multi;
            
            out = typecast(this.read(n),type);  
            
            if ~isempty(sz) && ~isempty(out)
                out = reshape(out,sz);
            end

        end
        
        function bytes = read( this, n )
            
            bytes = [];
            
            try
                if this.isAvailable(n)
                    bytes = org.apache.commons.io.IOUtils.toByteArray(this.inputStream,n);
                    this.BytesRemaining_ = this.BytesRemaining_ - length(bytes);
                else
                    this.throwOverRequesrError(n);                   
                end
            catch ME
                if isa(ME,'matlab.exception.JavaException')                    
                    errorObj = ME.ExceptionObject;                    
                    if isa(errorObj,'java.io.IOException')                        
                        this.throwOverRequesrError(n);
                    else                        
                        throwAsCaller(ME);
                    end
                    
                else
                    throwAsCaller(ME);
                end                
            end
            
        end
        
        function ME = throwOverRequesrError( this, n )
            ME = MException('DataReader:overrequest',...
                'To many bytes requested. Bytes remaining %d.  Bytes reuested %d',...
                this.BytesRemaining_, n);
            throw(ME);            
        end        

        function [ n, sz, ME ] = validateDimension( this, n )

            sz = [];
            ME = [];

            if isnumeric(n)  
                if any(isnan(n))
                    ME = getError(1); 
                    return
                end
                if isscalar(n)
                    if ~isfinite(n)
                        n = this.BytesRemaining_;
                    end                   
                else
                    if isvector(n)                
                        if (numel(n) == 2)
                            if any(isfinite(n))
                                handleInfinitSize();
                            else
                                sz = n;
                                n = prod(n);
                            end
                        else
                            ME = getError(2);
                        end
                    else
                        ME = getError(2);              
                    end                    
                end
            else
                n = [];
                ME = getError(2);        
            end
            
            function handleInfinitSize( )
                
                bytesRemaining = this.BytesRemaining_;
                
                finiteMask = isfinite(n);
                
                noFixed = n(finiteMask);
                
                noFlex = (bytesRemaining - rem(bytesRemaining,noFixed))/noFixed;
                
                if (noFlex > 0)
                
                    n(~finiteMask) = noFlex;
                
                    sz = n;
                
                    n = prod(n);
                else
                    % noFixed will allways be greater than bytesRemaining
                    % so will error in read by design
                    n = noFixed; 
                end

            end

            function ME = getError( type )
                
                switch type
                    case 1
                        msg = 'NaN values for size parameter is not allowed';
                    case 2                        
                        msg = 'The dimension input must be a numeric scalar or a two element row vector';        
                end
                
                ME = MException('DataReader:badinput',msg);
            end
        end 
    end
            
end


function [input,fileLength,fileName] = loadFromFile()

    input = [];
    fileLength = 0;
    
    [fileName,path] = uigetfile('*.*','Select File');

    if ~isequal(fileName,0)
        fileName = fullfile(path,fileName);
        file = java.io.File(fileName);
        fileLength = java.nio.file.Files.size(file);
        input = org.apache.commons.io.FileUtils.openInputStream(file);
    else
        fileName = [];
    end 
 
end