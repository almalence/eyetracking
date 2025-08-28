function entries = getEntry( this, which, includeDirectories )
% getEntry  
    narginchk(1,3);
    
    if nargin < 3
        includeDirectories = false;
    else
        try
            includeDirectories = logical(includeDirectories);
        catch
           ME = MException('ZipFile:NotLogical','The allowDirectories argrument must be converitible to a boolean');
           throwAsCaller(ME);            
        end
    end    
    
    if nargin == 1
        entries = getAll();
        return;
    end

    if isa(which,'io.Entry')
        entries = which;
        return
    elseif isempty(which)
        entries = getAll();
        return;
    end
    
    which = controllib.internal.util.hString2Char(which);

    if ischar(which) 
        entries = processCharInput();
    elseif iscell(which)
        
    elseif isnumeric(which)
        entries = processNumericInput();
    end 

    if isempty(entries)
       ME = MException('ZipFile:FindEntry','Could not find an archive entry match');
       throwAsCaller(ME);
    end

    function entries = getAll()
        entries = this.Entries;
        if ~includeDirectories
            notDirectoryMask = ~[this.Entries(:).IsDirectory];
            entries = entries(notDirectoryMask);
        end                
    end

    function entries = processCharInput()
        if strcmp(which,'all')
            entries = getAll();
        else
            mlocs = strfind(which,':');
            if isempty(mlocs)
                mask = endsWith(this.Files_,which);
                type = 1;
            elseif numel(mlocs) == 1
                if mlocs < numel(which) * 0.5
                    which = which(mlocs+1:end);
                    mask = startsWith(this.Files_,which);
                    type = 2;
                else
                    which = which(1:mlocs-1);
                    mask = endsWith(this.Files_,which); 
                    type = 1;
                end
            elseif numel(mlocs) == 2
                which = which(mlocs(1)+1:mlocs(2)-1);
                start = regexp(this.Files_,which,'start');
                mask = ~cellfun('isempty',start);
                type = 3;
            end

            if sum(mask) > 0
                entries = this.Entries(this.FilesIndexMap(mask));
            else
                switch type
                    case 1
                        ME = MException('ZipFile:FindEntry','Cannot locate any archive entries using endsWith and the pattern, ''%s''',which);
                    case 2
                        ME = MException('ZipFile:FindEntry','Cannot locate any archive entries using startsWith and the pattern, ''%s''',which);
                    case 3
                        ME = MException('ZipFile:FindEntry','Cannot locate any archive entries using the regular expression, ''%s''',which);
                end
                throwAsCaller(ME);
            end

        end        
    end

    function entries = processNumericInput()
        if isvector(which)

            minIdx = min(which);
            maxIdx = max(which);

            if (minIdx > 0) && (maxIdx <= numel(this.Files))
                entries = this.Entries(this.FilesIndexMap(which));
            else
                if minIdx < 1
                    msg = 'One or more of the provided indices is less than 1';
                elseif maxIdx > numel(this.Files)
                    msg = sprintf('One or more of the provided indices is greater than the number of entries (%d)',numel(this.Files));
                end
               ME = MException('ZipFile:FindEntry',msg);
               throwAsCaller(ME);
            end                  
        else
           ME = MException('ZipFile:FindEntry','Indices into the entry list must be a scalar or a vector');
           throwAsCaller(ME);
        end        
    end
        
end

