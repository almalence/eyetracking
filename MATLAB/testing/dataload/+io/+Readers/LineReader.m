classdef LineReader < io.Readers.AbstractReader
% LineReader - Object to read a compressed text file line by line
%
%   Syntax:

    properties( Dependent = true )
       
        LineNumber (1,1) double
    end

    
    methods
        function this = LineReader( entry, inputStream, buffSize )
            
            this@io.Readers.AbstractReader( entry, inputStream, buffSize );

        end
        
         function lineNumber = get.LineNumber( this )
             
            if this.AtEndOfFile_
                lineNumber = -1;
                return
            end
            
            try
                lineNumber = this.Reader.getLineNumber();
            catch
                lineNumber = -1;
            end
         end        

    end 
    
    methods
       
        function line = readLine( this )
            try                
                line = this.readLineImpl();
            catch ME
                io.Util.handleIOExceptions( ME );
            end
        end
        
        function lines = readLines( this, n )
            
            getAll = (nargin == 1);
            
            if ~getAll
                validateInteger( n, 'readLines' );

                lines = cell(n,1);

                for itr = 1:n

                    [line,isValid] = this.readLineImpl();
                    
                    if isValid
                        lines{itr} = line;
                    else
                        break
                    end

                end

                if ~isValid
                    lines(itr-1:end) = [];
                end
            else
                lines = [];
                
                isValid = true;

                while isValid
                    
                    [line,isValid] = this.readLineImpl();
                    
                    if isValid
                        
                        lines{end+1,1} = line; %#ok<AGROW>
                    end                    
                end
            end
            
        end

        function linesSkipped = skipLines( this, linesToSkip )
            
            this.validateInteger( linesToSkip, 'skipLines' );
            
            linesSkipped = linesToSkip;
            
            for itr = 1:linesToSkip
               
                [~,isValid] = this.readLineImpl();
                
                if ~isValid
                    linesSkipped = itr;
                    break
                end
            end
        end
    end

    methods(Access = protected)
        
        function createReader( this, inputStream, buffSize )
            
            this.Reader = java.io.LineNumberReader(java.io.InputStreamReader(inputStream),buffSize);
            
            this.Reader.setLineNumber(1);            
        end

        function [line,isValid] = readLineImpl( this )
            
            line = this.Reader.readLine();
            this.AtEndOfFile_ = isempty(line);
            isValid = ~this.AtEndOfFile_;
            if ~this.AtEndOfFile_
                line = char(line);
            end
            this.HasBeenRead_ = true;            
        end
        
    end
end


