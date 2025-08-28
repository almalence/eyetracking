classdef GeneralPurposeBit
    
    properties(Dependent = true)
        useEncyrption
        useDataDescriptor
        useStrongEncyrption
        useUTF8ForNames
        isLocalHeaderMasked
        slidingDictionarySize;
        numberOfShannonFanoTrees;
    end
    
    properties(Access = private)

        BITS = zeros([16,1],'uint8');
    end
    
    properties(Constant = true)
        
    end
    
    methods
        function this = GeneralPurposeBit( gpb )
            
            if nargin
                
                this.BITS = gpb;
                
            end
        end
        
        function boolean = get.useEncyrption( this )
            boolean = this.BITS(1);
        end
        
        function boolean = get.useDataDescriptor( this )
            boolean = this.BITS(4);
        end
        
        function boolean = get.useStrongEncyrption( this )
            boolean = this.BITS(6);
        end  
        
        function boolean = get.useUTF8ForNames( this )
            boolean = this.BITS(12);
        end 
        
        function boolean = get.isLocalHeaderMasked( this )
            boolean = this.BITS(14);
        end
        
        function sz = get.slidingDictionarySize( this )
            if this.BITS(2)
                sz = uint16(8192);
            else
                sz = uint16(4096);
            end
        end        
        
        function noTrees = get.numberOfShannonFanoTrees( this )
            if this.BITS(3)
                noTrees = uint8(3);
            else
                noTrees = uint8(2);
            end                
        end
    end
    
    methods(Access = public)
       
        function tf = isExtractable( this )
            tf = not(this.BITS(1) | this.BITS(6) | this.BITS(14));
        end
    end
end

