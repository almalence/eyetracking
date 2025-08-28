classdef ZipFile <  handle
% ZipFile Object that allows for individual extraction of archve entries from zip archive files.
%
%   Syntax:
%       zipFile = io.ZipFile(zipFileName) creates an io.ZipFile object from
%       the zip file, zipFileName.
%
%       zipFile = io.ZipFile()  creates an empty io.ZipFile object.
%       Setting the ArchiveFileName property to a valid zip file will load
%       the io.ZipFile object.
%
%   Properties:
%       ArchiveFileName     - Absoulte path to the zip file.
%       Files               - Cell array of archive entries
%       Directories         - Cell array of archives entries that are directories
%
%   Methods:
%       extract             - Extract an entry to file.
%       getDataReaderFor    - Returns a input stream reader for reading binary
%                             data directly to Matlab.
%       getLineReaderFor    - Returns an input stream reader for reading line
%                             terminated text directly to Matlab.
%       getInputStreamFor   - Returns an input stream for an entry.
%       getURLFor           - Returns the URL for an entry.
%       getAllImagesAsIcons - Returns a Java Array of ImageIcons of all
%                             recognized images in the archive file.
%       getComments         - Returns the archive file comments.
%       getSupportedModes   - Returns a array of supported decompression
%                             methods.
%       getEntry            - Returns a io.Entry object.
%
%   Notes:
%   - 

    
    properties(Dependent)
        % ArchiveFileName - Absoulte archive file name as returned by fopen
        ArchiveFileName     (1,:) char;     
        % Files - Cell array archive entry names that do not end with '/'
        Files               (:,1) cell; 
        % Directories - Cell array of directories
        Directories         (:,1) cell;
    end
    
    properties(Access = protected)
        % Private archive file name, files and directories
        ArchiveFile_         (1,:) char =  '';        
        Files_               (:,1) cell;        
        Directories_         (:,1) cell;
        % FilesIndexMap - Array length of Files_, which maps file names to
        % Entries indices.
        FilesIndexMap        (:,1) double;
        % Array of io.Entry objects
        Entries              (:,1) io.Entry;
        % Flag indicating wether to recreate the JavaZipFile object
        resetJavaZipFile     (1,1) logical = true;      
        % Raw bytes (uint8) from digital signature
        RawDigitalSignature_ (:,1) uint8;        
        % Array of supported decompress methods.
        SupportedModes_      (1,:) uint16;        
        % Archive comments
        ArchiveComments_     (1,:) char = ''; 
        % Flag indicating wether the io.Entry contructor should read the
        % local file header data
        ReadLocalFileHeaders_(1,1) logical = false;
        % A commons.compress ZipFile object. Created the first time that an
        % input stream is requested.
        JavaZipFile          = [];  
        % A com.mathworks StreamCopier object
        StreamCopier         = [];        
        % Raw bytes of the Zip64 extensible data
        Zip64ExtensibleData_ = [];   
        % Current io.ZipFile version
        Version = '1.0.0';
    end
    
    methods
        %==================================================================
        % Constructor
        %==================================================================
        function this = ZipFile( file )
           
            if nargin
                try
                    this.load(file);
                catch ME
                    throw(ME);
                end                
            else
                this.createUI();
                %TODO: UI?
            end
                
        end
        %==================================================================
        % Destructor
        %==================================================================        
        function delete( this )
            
           try this.JavaZipFile.close; catch; end 
        end
    end
    
    methods
        %==================================================================
        % Getters and setters
        %==================================================================         
        function set.ArchiveFileName( this, file )
            file = io.Util.validateZipFile(file);
            if strcmp(file,this.ArchiveFile_); return; end
            this.resetJavaZipFile = true;
            this.load(file);
        end
        
        function file = get.ArchiveFileName( this )
            file = this.ArchiveFile_;
        end
        
        function content = get.Files( this )
            content = this.Files_;
        end
        
        function directories = get.Directories( this )
            directories = this.Directories_;
        end 
    end
    
    methods(Access = public)
        %==================================================================
        % extract - Extract to file
        %==================================================================         
        function this = extract( this, whichEntries, outputFolder )
        % extract - Extracts the entry or entries defined by the argrument 
        % 'whichEntries' (see getEntry) to the folder 'outputFolder' or pwd.
        
            narginchk(1,3);
            % Use pwd for output directory unless a output directory is
            % provided, if outputFolder argrument is empty ([]), open folder dialog
            if nargin == 3
                if isempty(outputFolder)
                    outputFolder = uigetdir(pwd,'Select output folder');
                    
                    if isequal(outputFolder,0)
                        return
                    end
                end
                
                if ~exist(outputFolder,'dir')
                    mkdir(outputFolder)
                end                
            else
                
                outputFolder = pwd;
                if nargin == 1
                    whichEntries = [];
                end
            end

            entries = this.getEntry(whichEntries);
            
            for itr = 1:numel(entries)

                entry = entries(itr);

                [inputStream,msg] = this.getInputStreamForImpl( entry );

                if isempty(inputStream)
                    % TODO: should just warn
                    error('ZipFile:BadInputStream',msg);
                end                    

                outputFileName = fullfile(outputFolder,entry.FileName);

                jOutFile = java.io.File(outputFileName);

                outputStream = java.io.FileOutputStream(jOutFile);

                this.copyStream(inputStream, outputStream);
            end
            
        end
        %==================================================================
        % getDataReaderFor - Get a data reader for an entry
        %==================================================================         
        function dataReader = getDataReaderFor( this, whichEntry )
            
            entry = this.getEntry(whichEntry,false);
            entry = entry(1);

            dataReader = this.getReader( entry, 'data' );           
        end
        %==================================================================
        % getLineReaderFor - Get a line reader for an entry
        %==================================================================                 
        function lineReader = getLineReaderFor( this, whichEntry )
 
            entry = this.getEntry(whichEntry,false);
            entry = entry(1);
            
            lineReader = this.getReader( entry, 'line' );
        end
        %==================================================================
        % getTokenReaderFor - Get a token reader for an entry
        %==================================================================                 
        function lineReader = getTokenReaderFor( this, whichEntry )
 
            entry = this.getEntry(whichEntry,false);
            entry = entry(1);
            
            lineReader = this.getReader( entry, 'token' );
        end        
        %==================================================================
        % getInputStreamFor - Get a input stream for an entry
        %================================================================== 
        function inputStream = getInputStreamFor( this, whichEntry )
            
            entry = this.getEntry(whichEntry);
            entry = entry(1);
            
            inputStream = this.getInputStreamForImpl(entry);
        end
        %==================================================================
        % getURLFor - Get the URL(s) for an entry(s)
        %==================================================================         
        function URL = getURLFor( this, whichEntry )
            
            entries = this.getEntry(whichEntry,false);
            
            fullFileNames = arrayfun(@(x) x.getFullFileName, entries, 'UniformOutput',false)';
            
            URL = cellfun(@(x) io.Util.constructURL(x,this.ArchiveFile_),fullFileNames,'UniformOutput',false);
            
            if numel(URL) == 1
                URL = URL{1};
            end
            
        end
        %==================================================================
        % getAllImagesAsIcons - Get all images as ImageIcons
        %==================================================================          
        function icons = getAllImagesAsIcons( this )
            
            icons = [];
            
            files = this.Files;
            validExtensions = io.Util.getImageIOFileExtensions();
            
            idxs = find(contains(files,validExtensions));
            
            if ~isempty(idxs)
                entries = this.Entries(this.FilesIndexMap(idxs));
                URLs = this.getURLFor(entries);
                noURLs = numel(URLs);
                icons = javaArray('javax.swing.ImageIcon',noURLs);
               
                for itr = 1:noURLs
                    icons(itr) = javax.swing.ImageIcon(URLs{itr});
                end
            end
            
        end
        %==================================================================
        % getImage - Read a supported image file into Matlab directly
        %==================================================================          
        function [X,mapOrAlpha] = readImage( this, whichEntry, bgColor )
            
            supplyBG = (nargin == 3);
            
            entry = this.getEntry(whichEntry,false);
            entry = entry(1);
            
%             validExtensions = io.Util.getImageIOFileExtensions();
            
            %TODO: validate entry against valid extensions
            
            URL = io.Util.constructURL(entry.getFullFileName,this.ArchiveFile_);
            
            bi = javax.imageio.ImageIO.read(java.net.URL(URL));
            
            if supplyBG
                [X,mapOrAlpha] = io.internal.buffered2im(bi,bgColor);
            else
                [X,mapOrAlpha] = io.internal.buffered2im(bi);
            end
            
        end        
        %==================================================================
        % getDigitalSignature - Get the raw digital signature
        %==================================================================         
        function digSig = getDigitalSignature( this )
            digSig = this.RawDigitalSignature_;
        end
        %==================================================================
        % getSupportedModes - Get an array of supported modes
        %==================================================================         
        function modes = getSupportedModes( this )
            modes = this.SupportedModes_;
        end
        %==================================================================
        % getComments - Get file comments
        %==================================================================         
        function comments = getComments( this )
            if isempty(this.ArchiveComments_)
                comments = '';
            else
                comments = this.ArchiveComments_;
            end
        end
        %==================================================================
        % getEntry - Get an entry or entries (io.Entry)
        %==================================================================         
        entries = getEntry( this, which, allowDirectories ) 
    end

    methods( Access = protected )
        %==================================================================
        % load - Retrive archive entries
        %==================================================================         
        function load( this, file )

            try
                [entries,...
                 extensibleData,...
                 rawDigitalSig,...
                 archiveFile,...
                 archiveComments] = io.internal.getZipFileContent(file,this.ReadLocalFileHeaders_);

            catch ME
                throwAsCaller(ME);
            end
            
            this.Entries = entries;
            this.Zip64ExtensibleData_ = extensibleData;
            this.RawDigitalSignature_ = rawDigitalSig;
            this.ArchiveFile_ = archiveFile;
            this.ArchiveComments_ = archiveComments;
            
            notDirectoryMask = ~[entries(:).IsDirectory];
            
            this.Files_ = {entries(notDirectoryMask).FileName};
            
            this.Directories_ = {entries(~notDirectoryMask).Path};
            
            this.FilesIndexMap = find(notDirectoryMask)';  
            
            this.SupportedModes_ = io.Util.getSupportedModes;
        end
        %==================================================================
        % copyStream - Copies an inputstream to an outputstream
        %================================================================== 
        function copyStream( this, inputStream, outputStream )
           
            if isempty(this.StreamCopier)
                this.StreamCopier = com.mathworks.mlwidgets.io.InterruptibleStreamCopier.getInterruptibleStreamCopier; %#ok<*JAPIMATHWORKS>
            end
            
            try
                this.StreamCopier.copyStream(inputStream,outputStream);
                outputStream.close
            catch ME
                outputStream.close
                io.Util.handleIOExceptions(ME);
            end

        end
        %==================================================================
        % getInputStreamForImpl - getInputStream implementation
        %==================================================================         
        function [ inputStream , msg ] = getInputStreamForImpl( this, entry )
            
            msg = true;
                       
            if isempty(this.JavaZipFile) || this.resetJavaZipFile
                try
                    this.JavaZipFile =...
                        org.apache.commons.compress.archivers.zip.ZipFile(java.io.File(this.ArchiveFile_),[]);
                catch ME
                    io.Util.handleIOExceptions(ME);
                end
            end
            
            this.resetJavaZipFile = false;
            
            jEntry = this.JavaZipFile.getEntry(entry.getFullFileName);

            if isempty(jEntry)
                inputStream = [];
                msg = sprintf('Could not find a corresponding ZipArchive entry for %s',entry.FileName);
            elseif ~this.JavaZipFile.canReadEntryData(jEntry)
                inputStream = [];
                msg = sprintf('Found unsupported compression method for %s.\nTry updating org.apache.commons.compress to the latest version',entry.FileName);
            else
                try
                    inputStream = this.JavaZipFile.getInputStream(jEntry);
                catch ME
                    io.Util.handleIOExceptions(ME);
                end
            end           
        end
        %==================================================================
        % getReader - Gets a reader by type
        %==================================================================        
        function reader = getReader( this, entry, whichReader )

            [inputStream,msg] = this.getInputStreamForImpl(entry);
            
            if ~isempty(inputStream)
                switch whichReader
                    case 'data'
                        reader = io.Readers.DataReader( inputStream, entry, 4096);
                    case 'line'
                        reader = io.Readers.LineReader( entry, inputStream, 4096);
                    case 'token'
                        reader = io.Readers.TokenReader( entry, inputStream, 4096);
                end                        
            else
                ME = MException('ZipFile:BadInputStream',msg);
                throwAsCaller(ME);
            end           
        end
  
    end

end

function copyStream( inputStream, outputStream )

    try
        org.apache.commons.io.IOUtils.copy(inputStream,outputStream);
        outputStream.close();
    catch ME
        outputStream.close();
        throwAsCaller(ME);                
    end
end

function [panelfile,labelfile,browsebutton,panelselected,labelselected,...
    panelzip,labelzip] = getComponents(startPath,zipstr,selstr)

    border = javax.swing.UIManager.getLookAndFeel().getDefaults().getBorder('TextField.border');
    color = border.getLineColor;
    thick = border.getThickness;

    innerBorder = javax.swing.BorderFactory.createMatteBorder(thick,thick,0,thick,color);

    cls = 'com.mathworks.mwswing.MJPanel';
    panelfile = javaObjectEDT(cls);
    layout = com.jidesoft.swing.JideBoxLayout(panelfile,com.jidesoft.swing.JideBoxLayout.X_AXIS,2);
    panelfile.setLayout(layout);
    panelfile.setBorder(javax.swing.BorderFactory.createEmptyBorder(2,0,4,0));
    
    panelzip = javaObjectEDT(cls,java.awt.BorderLayout);

    panelselected = javaObjectEDT(cls,java.awt.BorderLayout);

    labelfile = createLabel(startPath,border);

    labelselected = createLabel(selstr,innerBorder);

    labelzip = createLabel(zipstr,innerBorder);
    
    cls = 'com.mathworks.mwswing.MJButton';
    icon = com.mathworks.common.icons.IconEnumerationUtils.getIcon('open_ts_16.png');
    browsebutton = handle(javaObjectEDT(cls,icon),'CallbackProperties');
    browsebutton.setName('browse');
    
    panelfile.add(labelfile,com.jidesoft.swing.JideBoxLayout.VARY);
    panelfile.add(browsebutton,com.jidesoft.swing.JideBoxLayout.FIX);
    
    panelselected.add(labelselected,'Center');
    
    panelzip.add(labelzip,'Center')
    
    function label = createLabel( str, border )
        
        cls = 'com.mathworks.mwswing.MJLabel';
        label = javaObjectEDT(cls,[' ',str]);
        label.setBorder(border);
        label.setMinimumSize(java.awt.Dimension(100,24));
    end

end