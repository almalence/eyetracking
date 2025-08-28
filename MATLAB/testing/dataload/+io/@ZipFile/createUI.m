function createUI( this )

    components = java.util.Hashtable();

    border = javax.swing.UIManager.getLookAndFeel().getDefaults().getBorder('TextField.border');

%     label = createLabel('',true,'zipfilenamelbl');

    button = createButton([],'open_ts_16.png',{@onSelectZipFile,this,components},'filebut');
    
    panelFileName = createPanel(true,[2,0,4,0]);
    cls = 'com.jidesoft.swing.LabeledTextField';
    ltf = handle(javaObjectEDT(cls),'CallbackProperties');
    ltf.setHintText('select zip file');
    ltf.setLabelText('Zip File : ');
    ltf.getTextField.setEnabled(false);
    components.put('zipfilenamelbl',ltf)
    
    addToPanel(panelFileName,ltf,'vary')
%     addToPanel(panelFileName,createLabel('Zip File', false),'fix')
%     addToPanel(panelFileName,label,'vary')
    addToPanel(panelFileName,button,'fix')

%     label = createLabel([pwd,filesep],true,'extractto');

    button = createButton([],'folders_ts_16.png',{@onSelectExtractTo,this,components},'folderbut');
    
    panelExtractTo = createPanel(true,[2,0,4,0]);
    cls = 'com.jidesoft.swing.LabeledTextField';
    ltf = handle(javaObjectEDT(cls),'CallbackProperties');
    ltf.setHintText('folder to extract to');
    ltf.setLabelText('Extract To : ');
    components.put('extractto',ltf)
    addToPanel(panelExtractTo,ltf,'vary')
%     addToPanel(panelExtractTo,createLabel('Extract To', false),'fix')
%     addToPanel(panelExtractTo,label,'vary')
    addToPanel(panelExtractTo,button,'fix')  
    
    panelTop = createPanel(false,[0,0,0,0]);    
    panelTop.add(panelFileName);
    panelTop.add(panelExtractTo);
    
    cls = 'com.jidesoft.list.DualList';    
    list = handle(javaObjectEDT(cls),'CallbackProperties');
    list.setRightButtonPanelVisible(false);
    list.setSelectionMode(com.jidesoft.list.DefaultDualListModel.REMOVE_SELECTION); 
    list.getOriginalListPane.getComponent(0).setBorder(javax.swing.BorderFactory.createEmptyBorder);
    list.getSelectedListPane.getComponent(0).setBorder(javax.swing.BorderFactory.createEmptyBorder);
    list.getOriginalListPane.setBorder(border);
    list.getSelectedListPane.setBorder(border); 
    list.getSelectedListPane.setBackground(java.awt.Color.WHITE);
    list.setBorder(javax.swing.BorderFactory.createEmptyBorder(4,0,2,0));
    components.put('duallist',list);

    label = createLabel('Zip File Content',true,'orglbl');
    label.setBorder(javax.swing.BorderFactory.createEmptyBorder);

    cls = 'com.jidesoft.list.QuickListFilterField';
    quickListField = handle(javaObjectEDT(cls),'CallbackProperties');
    quickListField.setPreferredSize(java.awt.Dimension(120,20))
    quickListField.setWildcardEnabled(true);
    quickListField.setBorder(javax.swing.BorderFactory.createEmptyBorder);    
    components.put('qlf',quickListField);
    
    panel = createPanel(true,[1,0,1,0]); 
    panel.setBackground(java.awt.Color.WHITE)
    addToPanel(panel,label,'vary')
    addToPanel(panel,quickListField,'fix')
    
    list.getOriginalListPane.add(panel,'North');

    label = createLabel('Files to be Extracted',true,'sellbl');  
    label.setPreferredSize(java.awt.Dimension(60,24))
    label.setBorder(javax.swing.BorderFactory.createEmptyBorder);

    list.getSelectedListPane.add(label,'North');

    button = createButton('Extract','import_ts_16.png',{@onProcess,this,components},'processbut');
    button.setEnabled(false);

    cls = 'com.mathworks.mwswing.MJProgressBar';    
%     cls = 'com.mathworks.toolbox.shared.computils.progress.widgets.CompUtilsJProgressBar';
    progressBar = handle(javaObjectEDT(cls),'CallbackProperties');
    progressBar.setBorderPainted(false);    
    components.put('pb',progressBar);
    
    panelBottom = createPanel(true,[8,0,2,0]);    
    addToPanel(panelBottom ,progressBar,'vary')
    addToPanel(panelBottom ,button,'fix')
    
    panelContent = javaObjectEDT(cls,java.awt.BorderLayout);
    panelContent.setBorder(javax.swing.BorderFactory.createEmptyBorder(8,8,8,8));    

    panelContent.add(panelTop,'North');
    panelContent.add(list,'Center');
    panelContent.add(panelBottom,'South');

    hFig = figure(...
        'Name','Zip File Import Tool',...
        'Tag','ZipFileImpotTool',...
        'MenuBar','none',...
        'Toolbar','none',...
        'DockControls','off',...
        'NumberTitle','off',...
        'Resize','on');
    
    hPanel = uipanel(...
        'Parent',hFig,...
        'BorderType','none',...
        'Units','norm',...
        'Position',[0,0,1,1]);    

    hgjavacomponent(...
        'Parent',hPanel,...
        'JavaPeer',panelContent,...
        'Units','normalized',...
        'Position',[0,0,1,1]); 
    
    function label = createLabel( str, hasBorder, name )
        
        if nargin < 3
            name = [];
        end
        
        cls = 'com.mathworks.widgets.SyntaxTextLabel';
%         cls = 'com.mathworks.mwswing.MJLabel';
        label = javaObjectEDT(cls,[' ',str],-1);
        if hasBorder
            label.setBorder(border);
        end
        label.setMinimumSize(java.awt.Dimension(100,24));
        
        if ~isempty(name)
            components.put(name,label);
        end
    end 

    function panel = createPanel( isHoriz, m )
                
        cls = 'com.mathworks.mwswing.MJPanel';
        panel = javaObjectEDT(cls);
        
        if isHoriz
            alignment = com.jidesoft.swing.JideBoxLayout.X_AXIS;
        else
            alignment = com.jidesoft.swing.JideBoxLayout.Y_AXIS;
        end
        layout = com.jidesoft.swing.JideBoxLayout(panel,alignment,2);
        panel.setLayout(layout);
        margin = javax.swing.BorderFactory.createEmptyBorder(m(1),m(2),m(3),m(4));
        panel.setBorder(margin);        
    end

    function button = createButton( text, icon, callback, name )

        cls = 'com.mathworks.mwswing.MJButton';
        icon = com.mathworks.common.icons.IconEnumerationUtils.getIcon(icon);
        if isempty(text)
            button = handle(javaObjectEDT(cls,icon),'CallbackProperties');
        else
            button = handle(javaObjectEDT(cls,text,icon),'CallbackProperties');
        end
        button.ActionPerformedCallback = callback; 
        components.put(name,button);
    end

    function addToPanel( panel, component, state )
        
        switch state
            case 'vary'
                state = com.jidesoft.swing.JideBoxLayout.VARY;
            case 'fix'
                state = com.jidesoft.swing.JideBoxLayout.FIX;
        end
        panel.add(component,state);
    end
end

function onSelectZipFile( src, evnt, this, components )

    path = getCurrentPath(components);

    zipFilePath = io.Util.uigetzipfile(path);

    if isempty(zipFilePath); return; end

    this.load(zipFilePath);

    [~,zipFileName,ext] = fileparts(zipFilePath);
    zipFileName = [zipFileName,ext];
    files = this.Files;
    [~,~,e] = fileparts(files);
    [~,I] = sort(e);            
    files = files(I);
    
    components.get('zipfilenamelbl').setText([' ',zipFileName]);
    
    cls = 'com.jidesoft.list.DefaultDualListModel';
    model = javaObjectEDT(cls);

    for itr = 1:numel(files)
        model.addElement(sprintf(' %d. %s',I(itr),files{itr}));
    end 
    
    model = handle(model,'CallbackProperties');
    model.ValueChangedCallback = {@onValueChanged,components};
    
    components.get('qlf').setListModel(model);

    components.get('duallist').setModel(model);
    components.get('duallist').getOriginalList.setModel(components.get('qlf').getDisplayListModel);
    
    components.get('orglbl').setText(sprintf(' Zip File Content (%d)',numel(files)));
    
end

function onProcess( src, evnt, this, components )

    names = components.get('duallist').getSelectedValues();
    
    path = getCurrentPath(components);
    
    noNames = numel(names);
    
    doPB = false;
    
    if noNames > 1
        progressBar = components.get('pb');
        progressBar.setMinimum(0);
        progressBar.setMaximum(noNames);
        progressBar.setValue(0);
        progressBar.setStringPainted(true);
        doPB = true;
    end
    
    for itr = 1:noNames
        name = char(names(itr));
        name = strtrim(name(find(name == 46,1,'first')+2:end));
        if doPB
            progressBar.setValue(itr);
            progressBar.setString(sprintf('Extracting... %d/%d',itr,noNames));
        end
        this.extract(name,path);
    end   
    
    if doPB
        progressBar.setStringPainted(false);
        progressBar.setValue(0);
    end

end

function onValueChanged( src, evnt, components )

    components.get('processbut').setEnabled(~src.isSelectionEmpty);
    
    selCnt = numel(components.get('duallist').getSelectedIndices);
    orgCnt = numel(components.get('duallist').getUnselectedIndices);
    
    components.get('orglbl').setText(sprintf(' Zip File Content (%d)',orgCnt));
    if selCnt > 0
        components.get('sellbl').setText(sprintf(' Files to be Extracted (%d)',selCnt));
    else
        components.get('sellbl').setText(' Files to be Extracted');
    end
end

function onSelectExtractTo( src, evnt, this, components )

    path = getCurrentPath(components);
    
    selpath = uigetdir(path,'Select directory to extract files to');
    
    if ~isequal(selpath,0)
        components.get('extractto').setText([' ',selpath,filesep]);
    end
end

function path = getCurrentPath( components )

    path = strtrim(char(components.get('extractto').getText));
end