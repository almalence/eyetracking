classdef AbstractReader < handle

    
    properties( Dependent = true )
        FileName                (1,:) char
        
        Entry                   (1,1) io.Entry
        
        FileSize                (1,1) double
        
        CompressedSize          (1,1) doube
    end
    
    properties( Access = protected )
        
        Reader        
       
        Entry_                  (1,1) io.Entry;    
        
        AtEndOfFile_            (1,1) logical =  false;
        
        HasBeenRead_            (1,1) logical =  false;
        
    end
    
    methods
        function this = AbstractReader( entry, inputStream, buffSize )

            this.Entry_ = entry;
            
            uncompressedSize = entry.UncompressedSize;
            
            if nargin < 3
                buffSize = 4096;
            else
               this.validateInteger( buffSize, 'Constructor' );
            end
            
            if uncompressedSize < buffSize
                buffSize = uncompressedSize;
            end
            
            this.createReader( inputStream, buffSize )
            
        end
        
        function delete( this )
           
            try this.Reader.close(); catch; end
        end

    end
    
    methods
       
        function fileName = get.FileName( this )
           fileName = this.Entry_.FileName; 
        end
        
        function entry = get.Entry( this )
            entry = this.Entry_; 
        end 
         
        function fileSize = get.FileSize( this )
            fileSize = this.Entry_.UncompressedSize; 
        end
        
        function compressedSize = get.CompressedSize( this )
            compressedSize = this.Entry_.CompressedSize; 
        end        
        
        function reader = getReader( this )
            reader = this.Reader;
        end
    end 
    
    methods(Access = protected)
        
        function createReader( this, inputStream, buffSize ) %#ok<*INUSD>
            
        end
    end
  
    methods(Access = protected, Static = true)
        function validateInteger( n, fcnName )

            try
                validateattributes(n,{'numeric'},{'scalar','integer'},fcnName)
            catch ME
                throwAsCaller(ME)
            end
        end  
        
        function validateBoolean( n, fcnName )

            try
                validateattributes(n,{'nlogical'},{'scalar'},fcnName)
            catch ME
                throwAsCaller(ME)
            end
        end         
    end
end



