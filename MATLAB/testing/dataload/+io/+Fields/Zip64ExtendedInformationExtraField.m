classdef Zip64ExtendedInformationExtraField <  io.Fields.Field

    methods
        function this = Zip64ExtendedInformationExtraField( id, dataBytes  )
            this = this@io.Fields.Field(id,dataBytes);
            this.Name_ = 'Zip64ExtendedInformationExtraField';
        end
        
        function [msg,compressedSize,uncompressedSize,offset,diskNumber] =...
                update(this,compressedSize,uncompressedSize,offset,diskNumber)
            
            msg = [];
            neededFlags = [compressedSize,uncompressedSize,offset,diskNumber] < 0;
            
            if any(neededFlags)
%                requiredSize = ( neededFlags * io.Util.Z64_EXT_SIZE ) + 4;
                requiredSize = ( neededFlags * io.Util.Z64_EXT_SIZE );
                fieldSize = numel(this.RawData_);
                if fieldSize < requiredSize
                    msg = 'The Zip64 extended infomation field is not large enough to contain all required information';
                    return;
                end  
                newVals = zeros([1,4]);
                noBytesPerField = neededFlags .* io.Util.Z64_EXT_SIZE';
                ptr = 1;
                for itr = 1:4
                    if neededFlags(itr)
                        noBytes = noBytesPerField(itr);
                        if noBytes == 8
                            clazz = 'uint64';
                        else
                            clazz = 'uint32';
                        end
                            
                        newVals(itr) = double(typecast(this.RawData_(ptr:ptr+noBytes-1),clazz));
                        ptr = ptr + noBytes;
                    end
                end

                if isequal(compressedSize,-1)
                    compressedSize = newVals(1);
                else
                    compressedSize = double(compressedSize);
                end

                if isequal(uncompressedSize,-1)
                    uncompressedSize = newVals(2);
                else
                    uncompressedSize = double(uncompressedSize);
                end

                if isequal(offset,-1)
                    offset = newVals(3);
                else
                    offset = double(offset);
                end  

                if isequal(diskNumber,-1)
                    diskNumber = newVals(4);
                else
                    diskNumber = double(diskNumber);
                end                 
            else
                compressedSize = double(compressedSize);
                uncompressedSize = double(uncompressedSize);
                offset = double(offset);
                diskNumber = double(diskNumber);
            end
        end
    end
end

