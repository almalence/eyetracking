classdef Field < matlab.mixin.Heterogeneous
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Dependent = true)
       Id
       Name
    end
    properties(Access = protected)
        Id_(1,1) uint16
        RawData_(:,1) uint8
        Name_ = 'Generic'
    end
    
    methods
        function this = Field( id, dataBytes  )
            narginchk(0,2)
            if nargin
                this.Id_ = id;
                this.RawData_ = dataBytes;
            end
        end
        
        function value = get.Id( this )
            value = this.Id_;
        end
        
        function value = get.Name( this )
            value = this.Name_;
        end
    end
    
    
end

