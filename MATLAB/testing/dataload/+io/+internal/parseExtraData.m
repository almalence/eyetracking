function fields = parseExtraData( bytes, has4ByteDataSize )

    if nargin == 1
        has4ByteDataSize = false;
    end

    if ~has4ByteDataSize

        inc = 3;
        clazz = 'uint16';
    else
        inc = 7;
        clazz = 'uint32';
    end

    noBytes = numel(bytes);

    if noBytes < 4
        fields = [];
        return
    end

    doLoop = true;

    k = 1;
    ptr = 1;
    while doLoop

        id = typecast(bytes(ptr:ptr+1),'uint16');
        fsize = typecast(bytes(ptr+2:ptr+inc),clazz);
        data = bytes(ptr+inc+1:ptr+inc+fsize);
        switch id
            case 1
                fields(k) = io.Fields.Zip64ExtendedInformationExtraField(id,data); %#ok<AGROW>
            otherwise
                fields(k) = io.Fields.Field(id,data); %#ok<AGROW>
        end
        ptr = ptr + 4 + fsize;
        k = k + 1;
        doLoop = ptr < noBytes;
    end
end
