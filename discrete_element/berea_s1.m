function h1 = berea_s1()
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
% This problem is solved by saving the output as a FIG-file.
% 
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.
% 
% NOTE: certain newer features in MATLAB may not have been saved in this
% M-file due to limitations of this format, which has been superseded by
% FIG-files.  Figures which have been annotated using the plot editor tools
% are incompatible with the M-file/MAT-file format, and should be saved as
% FIG-files.


load berea_s1.mat


h1 = figure(...
'Color',[1 1 1],...
'Colormap',get(0,'defaultfigureColormap'),...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'Position',get(0,'defaultfigurePosition'));

appdata = [];
appdata.PlotColorIndex = mat{1};
appdata.PlotLineStyleIndex = mat{2};
appdata.PlotHoldStyle = mat{3};
appdata.LegendColorbarText = mat{4};
appdata.LegendColorbarInnerList = mat{5};
appdata.LegendColorbarOuterList = mat{6};
appdata.inLayout = mat{7};
appdata.LegendComputePosCache = mat{8};
appdata.LegendPeerHandle = mat{9};

h2 = axes(...
'Parent',h1,...
'Box','on',...
'CameraPosition',[0.0515463917525775 125000000 17.3205080756888],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'NextPlot','add',...
'XColor',get(0,'defaultaxesXColor'),...
'XLim',[0 0.103092783505155],...
'XLimMode','manual',...
'YColor',get(0,'defaultaxesYColor'),...
'ZColor',get(0,'defaultaxesZColor'),...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h3 = get(h2,'title');

set(h3,...
'Parent',h2,...
'Units','data',...
'FontUnits','points',...
'BackgroundColor','none',...
'Color',[0 0 0],...
'DisplayName',blanks(0),...
'EdgeColor','none',...
'EraseMode','normal',...
'DVIMode','auto',...
'FontAngle','normal',...
'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','normal',...
'HorizontalAlignment','center',...
'LineStyle','-',...
'LineWidth',0.5,...
'Margin',2,...
'Position',[0.0514276212646683 254751461.988304 1.00010919874411],...
'Rotation',0,...
'String',blanks(0),...
'Interpreter','tex',...
'VerticalAlignment','bottom',...
'ButtonDownFcn',[],...
'CreateFcn', {@local_CreateFcn, [], ''} ,...
'DeleteFcn',[],...
'BusyAction','queue',...
'HandleVisibility','off',...
'HelpTopicKey',blanks(0),...
'HitTest','on',...
'Interruptible','on',...
'SelectionHighlight','on',...
'Serializable','on',...
'Tag',blanks(0),...
'UserData',[],...
'Visible','on',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'IncludeRenderer','on',...
'Clipping','off');

h4 = get(h2,'xlabel');

set(h4,...
'Parent',h2,...
'Units','data',...
'FontUnits','points',...
'BackgroundColor','none',...
'Color',[0 0 0],...
'DisplayName',blanks(0),...
'EdgeColor','none',...
'EraseMode','normal',...
'DVIMode','auto',...
'FontAngle','normal',...
'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','normal',...
'HorizontalAlignment','center',...
'LineStyle','-',...
'LineWidth',0.5,...
'Margin',2,...
'Position',[0.0514276212646683 -17178362.5730994 1.00010919874411],...
'Rotation',0,...
'String','Axial Strain',...
'Interpreter','tex',...
'VerticalAlignment','cap',...
'ButtonDownFcn',[],...
'CreateFcn', {@local_CreateFcn, [], ''} ,...
'DeleteFcn',[],...
'BusyAction','queue',...
'HandleVisibility','off',...
'HelpTopicKey',blanks(0),...
'HitTest','on',...
'Interruptible','on',...
'SelectionHighlight','on',...
'Serializable','on',...
'Tag',blanks(0),...
'UserData',[],...
'Visible','on',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'IncludeRenderer','on',...
'Clipping','off');

h5 = get(h2,'ylabel');

set(h5,...
'Parent',h2,...
'Units','data',...
'FontUnits','points',...
'BackgroundColor','none',...
'Color',[0 0 0],...
'DisplayName',blanks(0),...
'EdgeColor','none',...
'EraseMode','normal',...
'DVIMode','auto',...
'FontAngle','normal',...
'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','normal',...
'HorizontalAlignment','center',...
'LineStyle','-',...
'LineWidth',0.5,...
'Margin',2,...
'Position',[-0.00676991781082238 123903508.77193 1.00010919874411],...
'Rotation',90,...
'String','Axial Stress (MPa)',...
'Interpreter','tex',...
'VerticalAlignment','bottom',...
'ButtonDownFcn',[],...
'CreateFcn', {@local_CreateFcn, [], ''} ,...
'DeleteFcn',[],...
'BusyAction','queue',...
'HandleVisibility','off',...
'HelpTopicKey',blanks(0),...
'HitTest','on',...
'Interruptible','on',...
'SelectionHighlight','on',...
'Serializable','on',...
'Tag',blanks(0),...
'UserData',[],...
'Visible','on',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'IncludeRenderer','on',...
'Clipping','off');

h6 = get(h2,'zlabel');

set(h6,...
'Parent',h2,...
'Units','data',...
'FontUnits','points',...
'BackgroundColor','none',...
'Color',[0 0 0],...
'DisplayName',blanks(0),...
'EdgeColor','none',...
'EraseMode','normal',...
'DVIMode','auto',...
'FontAngle','normal',...
'FontName','Helvetica',...
'FontSize',10,...
'FontWeight','normal',...
'HorizontalAlignment','right',...
'LineStyle','-',...
'LineWidth',0.5,...
'Margin',2,...
'Position',[-0.0174592617226472 272295321.637427 1.00010919874411],...
'Rotation',0,...
'String',blanks(0),...
'Interpreter','tex',...
'VerticalAlignment','middle',...
'ButtonDownFcn',[],...
'CreateFcn', {@local_CreateFcn, [], ''} ,...
'DeleteFcn',[],...
'BusyAction','queue',...
'HandleVisibility','off',...
'HelpTopicKey',blanks(0),...
'HitTest','on',...
'Interruptible','on',...
'SelectionHighlight','on',...
'Serializable','on',...
'Tag',blanks(0),...
'UserData',[],...
'Visible','off',...
'XLimInclude','on',...
'YLimInclude','on',...
'ZLimInclude','on',...
'CLimInclude','on',...
'ALimInclude','on',...
'IncludeRenderer','on',...
'Clipping','off');

appdata = [];
appdata.legend_texthandle = mat{10};
appdata.legend_linetokenhandle = mat{11};
appdata.legend_linemarkertokenhandle = mat{12};

h7 = graph2d.lineseries(...
'Parent',h2,...
'DisplayName','0 MPa',...
'Color',[0.368627450980392 0.309803921568627 0.635294117647059],...
'LineWidth',2,...
'XData',mat{13},...
'YData',mat{14},...
'XDataJitter',0,...
'XDataMode','manual',...
'ObeyXDataMode','auto',...
'XDataSource',blanks(0),...
'YDataSource',blanks(0),...
'ZDataSource',blanks(0),...
'CodeGenColorMode','manual',...
'CodeGenLineStyleMode','auto',...
'CodeGenMarkerMode','auto',...
'SwitchProps',{  'LineWidth'; 'LineStyle'; 'Color'; 'MarkerEdgeColor'; 'MarkerFaceColor'; 'Marker'; 'MarkerSize'; 'Visible'; 'XDataSource'; 'YDataSource'; 'ZDataSource'; 'DisplayName' },...
'OldSwitchProps',[],...
'OldSwitchVals',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.legend_texthandle = mat{15};
appdata.legend_linetokenhandle = mat{16};
appdata.legend_linemarkertokenhandle = mat{17};

h8 = graph2d.lineseries(...
'Parent',h2,...
'DisplayName','2 MPa',...
'Color',[0.200519399225414 0.55934234524082 0.738031469377875],...
'LineWidth',2,...
'XData',mat{18},...
'YData',mat{19},...
'XDataJitter',0,...
'XDataMode','manual',...
'ObeyXDataMode','auto',...
'XDataSource',blanks(0),...
'YDataSource',blanks(0),...
'ZDataSource',blanks(0),...
'CodeGenColorMode','manual',...
'CodeGenLineStyleMode','auto',...
'CodeGenMarkerMode','auto',...
'SwitchProps',{  'LineWidth'; 'LineStyle'; 'Color'; 'MarkerEdgeColor'; 'MarkerFaceColor'; 'Marker'; 'MarkerSize'; 'Visible'; 'XDataSource'; 'YDataSource'; 'ZDataSource'; 'DisplayName' },...
'OldSwitchProps',[],...
'OldSwitchVals',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.legend_texthandle = mat{20};
appdata.legend_linetokenhandle = mat{21};
appdata.legend_linemarkertokenhandle = mat{22};

h9 = graph2d.lineseries(...
'Parent',h2,...
'DisplayName','5 MPa',...
'Color',[0.455832079626395 0.789744115511128 0.645829808804391],...
'LineWidth',2,...
'XData',mat{23},...
'YData',mat{24},...
'XDataJitter',0,...
'XDataMode','manual',...
'ObeyXDataMode','auto',...
'XDataSource',blanks(0),...
'YDataSource',blanks(0),...
'ZDataSource',blanks(0),...
'CodeGenColorMode','manual',...
'CodeGenLineStyleMode','auto',...
'CodeGenMarkerMode','auto',...
'SwitchProps',{  'LineWidth'; 'LineStyle'; 'Color'; 'MarkerEdgeColor'; 'MarkerFaceColor'; 'Marker'; 'MarkerSize'; 'Visible'; 'XDataSource'; 'YDataSource'; 'ZDataSource'; 'DisplayName' },...
'OldSwitchProps',[],...
'OldSwitchVals',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.legend_texthandle = mat{25};
appdata.legend_linetokenhandle = mat{26};
appdata.legend_linemarkertokenhandle = mat{27};

h10 = graph2d.lineseries(...
'Parent',h2,...
'DisplayName','10 MPa',...
'Color',[0.761635738765942 0.905831090606177 0.629864253393665],...
'LineWidth',2,...
'XData',mat{28},...
'YData',mat{29},...
'XDataJitter',0,...
'XDataMode','manual',...
'ObeyXDataMode','auto',...
'XDataSource',blanks(0),...
'YDataSource',blanks(0),...
'ZDataSource',blanks(0),...
'CodeGenColorMode','manual',...
'CodeGenLineStyleMode','auto',...
'CodeGenMarkerMode','auto',...
'SwitchProps',{  'LineWidth'; 'LineStyle'; 'Color'; 'MarkerEdgeColor'; 'MarkerFaceColor'; 'Marker'; 'MarkerSize'; 'Visible'; 'XDataSource'; 'YDataSource'; 'ZDataSource'; 'DisplayName' },...
'OldSwitchProps',[],...
'OldSwitchVals',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.legend_texthandle = mat{30};
appdata.legend_linetokenhandle = mat{31};
appdata.legend_linemarkertokenhandle = mat{32};

h11 = graph2d.lineseries(...
'Parent',h2,...
'DisplayName','15 MPa',...
'Color',[0.92774517117992 0.958344124115124 0.644238952096614],...
'LineWidth',2,...
'XData',mat{33},...
'YData',mat{34},...
'XDataJitter',0,...
'XDataMode','manual',...
'ObeyXDataMode','auto',...
'XDataSource',blanks(0),...
'YDataSource',blanks(0),...
'ZDataSource',blanks(0),...
'CodeGenColorMode','manual',...
'CodeGenLineStyleMode','auto',...
'CodeGenMarkerMode','auto',...
'SwitchProps',{  'LineWidth'; 'LineStyle'; 'Color'; 'MarkerEdgeColor'; 'MarkerFaceColor'; 'Marker'; 'MarkerSize'; 'Visible'; 'XDataSource'; 'YDataSource'; 'ZDataSource'; 'DisplayName' },...
'OldSwitchProps',[],...
'OldSwitchVals',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.legend_texthandle = mat{35};
appdata.legend_linetokenhandle = mat{36};
appdata.legend_linemarkertokenhandle = mat{37};

h12 = graph2d.lineseries(...
'Parent',h2,...
'DisplayName','20 MPa',...
'Color',[0.982025247227377 0.920603629584506 0.637231663608876],...
'LineWidth',2,...
'XData',mat{38},...
'YData',mat{39},...
'XDataJitter',0,...
'XDataMode','manual',...
'ObeyXDataMode','auto',...
'XDataSource',blanks(0),...
'YDataSource',blanks(0),...
'ZDataSource',blanks(0),...
'CodeGenColorMode','manual',...
'CodeGenLineStyleMode','auto',...
'CodeGenMarkerMode','auto',...
'SwitchProps',{  'LineWidth'; 'LineStyle'; 'Color'; 'MarkerEdgeColor'; 'MarkerFaceColor'; 'Marker'; 'MarkerSize'; 'Visible'; 'XDataSource'; 'YDataSource'; 'ZDataSource'; 'DisplayName' },...
'OldSwitchProps',[],...
'OldSwitchVals',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.legend_texthandle = mat{40};
appdata.legend_linetokenhandle = mat{41};
appdata.legend_linemarkertokenhandle = mat{42};

h13 = graph2d.lineseries(...
'Parent',h2,...
'DisplayName','25 MPa',...
'Color',[0.994219317356572 0.758258112279141 0.431162294534286],...
'LineWidth',2,...
'XData',mat{43},...
'YData',mat{44},...
'XDataJitter',0,...
'XDataMode','manual',...
'ObeyXDataMode','auto',...
'XDataSource',blanks(0),...
'YDataSource',blanks(0),...
'ZDataSource',blanks(0),...
'CodeGenColorMode','manual',...
'CodeGenLineStyleMode','auto',...
'CodeGenMarkerMode','auto',...
'SwitchProps',{  'LineWidth'; 'LineStyle'; 'Color'; 'MarkerEdgeColor'; 'MarkerFaceColor'; 'Marker'; 'MarkerSize'; 'Visible'; 'XDataSource'; 'YDataSource'; 'ZDataSource'; 'DisplayName' },...
'OldSwitchProps',[],...
'OldSwitchVals',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.legend_texthandle = mat{45};
appdata.legend_linetokenhandle = mat{46};
appdata.legend_linemarkertokenhandle = mat{47};

h14 = graph2d.lineseries(...
'Parent',h2,...
'DisplayName','30 MPa',...
'Color',[0.968399903171145 0.479865737670985 0.272320395922429],...
'LineWidth',2,...
'XData',mat{48},...
'YData',mat{49},...
'XDataJitter',0,...
'XDataMode','manual',...
'ObeyXDataMode','auto',...
'XDataSource',blanks(0),...
'YDataSource',blanks(0),...
'ZDataSource',blanks(0),...
'CodeGenColorMode','manual',...
'CodeGenLineStyleMode','auto',...
'CodeGenMarkerMode','auto',...
'SwitchProps',{  'LineWidth'; 'LineStyle'; 'Color'; 'MarkerEdgeColor'; 'MarkerFaceColor'; 'Marker'; 'MarkerSize'; 'Visible'; 'XDataSource'; 'YDataSource'; 'ZDataSource'; 'DisplayName' },...
'OldSwitchProps',[],...
'OldSwitchVals',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.legend_texthandle = mat{50};
appdata.legend_linetokenhandle = mat{51};
appdata.legend_linemarkertokenhandle = mat{52};

h15 = graph2d.lineseries(...
'Parent',h2,...
'DisplayName','40 MPa',...
'Color',[0.852513916020359 0.2653890494876 0.308190107318648],...
'LineWidth',2,...
'XData',mat{53},...
'YData',mat{54},...
'XDataJitter',0,...
'XDataMode','manual',...
'ObeyXDataMode','auto',...
'XDataSource',blanks(0),...
'YDataSource',blanks(0),...
'ZDataSource',blanks(0),...
'CodeGenColorMode','manual',...
'CodeGenLineStyleMode','auto',...
'CodeGenMarkerMode','auto',...
'SwitchProps',{  'LineWidth'; 'LineStyle'; 'Color'; 'MarkerEdgeColor'; 'MarkerFaceColor'; 'Marker'; 'MarkerSize'; 'Visible'; 'XDataSource'; 'YDataSource'; 'ZDataSource'; 'DisplayName' },...
'OldSwitchProps',[],...
'OldSwitchVals',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.legend_texthandle = mat{55};
appdata.legend_linetokenhandle = mat{56};
appdata.legend_linemarkertokenhandle = mat{57};

h16 = graph2d.lineseries(...
'Parent',h2,...
'DisplayName','50 MPa',...
'Color',[0.619607843137255 0.00392156862745098 0.258823529411765],...
'LineWidth',2,...
'XData',mat{58},...
'YData',mat{59},...
'XDataJitter',0,...
'XDataMode','manual',...
'ObeyXDataMode','auto',...
'XDataSource',blanks(0),...
'YDataSource',blanks(0),...
'ZDataSource',blanks(0),...
'CodeGenColorMode','manual',...
'CodeGenLineStyleMode','auto',...
'CodeGenMarkerMode','auto',...
'SwitchProps',{  'LineWidth'; 'LineStyle'; 'Color'; 'MarkerEdgeColor'; 'MarkerFaceColor'; 'Marker'; 'MarkerSize'; 'Visible'; 'XDataSource'; 'YDataSource'; 'ZDataSource'; 'DisplayName' },...
'OldSwitchProps',[],...
'OldSwitchVals',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h17 = text(...
'Parent',h2,...
'Color',get(0,'defaulttextColor'),...
'HandleVisibility','off',...
'Tag','LegendDeleteProxy',...
'Visible','off');

h18 = text(...
'Parent',h2,...
'Units','normalized',...
'Color',get(0,'defaulttextColor'),...
'HorizontalAlignment','center',...
'String',' ',...
'DeleteFcn','legendcolorbarlayout(get(gcbo,''Parent''),''remove'')',...
'HandleVisibility','off',...
'HitTest','off',...
'Tag','LegendColorbarLayout');

h19 = text(...
'Parent',h2,...
'Units','normalized',...
'Color',get(0,'defaulttextColor'),...
'HorizontalAlignment','center',...
'Position',[1 1 0],...
'String',' ',...
'HandleVisibility','off',...
'HitTest','off',...
'Tag','LegendColorbarLayout');

appdata = [];
appdata.LegendTempText = mat{60};
appdata.MWBYPASS_title = mat{61};
appdata.MWBYPASS_xlabel = mat{62};
appdata.MWBYPASS_ylabel = mat{63};
appdata.MWBYPASS_zlabel = mat{64};
appdata.NonDataObject = mat{65};
appdata.PostDeserializeFcn = mat{66};
appdata.LegendOldSize = mat{67};

h20 = scribe.legend(...
'Parent',h1,...
'Position',[0.143154761904762 0.475 0.176785714285714 0.431746031746032],...
'Box','on',...
'CameraPosition',[0.5 0.5 17.3205080756888],...
'CameraPositionMode','auto',...
'CLim',[0 1],...
'CLimMode','manual',...
'Color',[1 1 1],...
'ColorOrder',mat{68},...
'DrawMode','fast',...
'FontSize',[],...
'LooseInset',[0 0 0 0],...
'NextPlot','add',...
'XColor',[0 0 0],...
'XLim',[0 1],...
'XLimMode','manual',...
'XTick',-1,...
'XTickLabel',blanks(0),...
'XTickLabelMode','manual',...
'XTickMode','manual',...
'YColor',[0 0 0],...
'YLim',[0 1],...
'YLimMode','manual',...
'YTick',-1,...
'YTickLabel',blanks(0),...
'YTickLabelMode','manual',...
'YTickMode','manual',...
'ZColor',[0 0 0],...
'ButtonDownFcn',mat{69},...
'Interruptible','off',...
'Tag','legend',...
'UserData',struct(...
    'PlotHandle', [], ...
    'legendpos', 2, ...
    'LegendPosition', [0.143154761904762 0.475 0.176785714285714 0.431746031746032], ...
    'LabelHandles', [], ...
    'handles', [], ...
    'lstrings', { {  '0 MPa'; '2 MPa'; '5 MPa'; '10 MPa'; '15 MPa'; '20 MPa'; '25 MPa'; '30 MPa'; '40 MPa'; '50 MPa' } }, ...
    'LegendHandle', []),...
'PlotChildListen','off',...
'Image',[],...
'Location','NorthWest',...
'Orientation','vertical',...
'ItemTextLocation','left',...
'ItemTokenSize',[30 18],...
'StyleLegend','on',...
'EdgeColor',[0 0 0],...
'ObserveText','on',...
'TextColor',[0 0 0],...
'Interpreter','tex',...
'String',{  '0 MPa' '2 MPa' '5 MPa' '10 MPa' '15 MPa' '20 MPa' '25 MPa' '30 MPa' '40 MPa' '50 MPa' },...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h21 = text(...
'Parent',h20,...
'Units','points',...
'Color',get(0,'defaulttextColor'),...
'FontSize',[],...
'Margin',0.01,...
'String','50 MPa',...
'HandleVisibility','off',...
'Tag','temphackytext',...
'Visible','off');

appdata = [];
appdata.Listeners = mat{70};

h22 = text(...
'Parent',h20,...
'Units','normalized',...
'Color',get(0,'defaulttextColor'),...
'FontSize',[],...
'Position',[0.525252525252525 0.938419117647059 0],...
'String','0 MPa',...
'ButtonDownFcn',mat{71},...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h23 = line(...
'Parent',h20,...
'Color',[0.368627450980392 0.309803921568627 0.635294117647059],...
'LineWidth',2,...
'XData',[0.0808080808080808 0.484848484848485],...
'YData',[0.938419117647059 0.938419117647059],...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'Tag','0 MPa');

h24 = line(...
'Parent',h20,...
'Color',[0.368627450980392 0.309803921568627 0.635294117647059],...
'LineStyle','none',...
'LineWidth',2,...
'XData',0.282828282828283,...
'YData',0.938419117647059,...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off');

appdata = [];
appdata.Listeners = mat{72};

h25 = text(...
'Parent',h20,...
'Units','normalized',...
'Color',get(0,'defaulttextColor'),...
'FontSize',[],...
'Position',[0.525252525252525 0.840992647058824 0],...
'String','2 MPa',...
'ButtonDownFcn',mat{73},...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h26 = line(...
'Parent',h20,...
'Color',[0.200519399225414 0.55934234524082 0.738031469377875],...
'LineWidth',2,...
'XData',[0.0808080808080808 0.484848484848485],...
'YData',[0.840992647058824 0.840992647058824],...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'Tag','2 MPa');

h27 = line(...
'Parent',h20,...
'Color',[0.200519399225414 0.55934234524082 0.738031469377875],...
'LineStyle','none',...
'LineWidth',2,...
'XData',0.282828282828283,...
'YData',0.840992647058824,...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off');

appdata = [];
appdata.Listeners = mat{74};

h28 = text(...
'Parent',h20,...
'Units','normalized',...
'Color',get(0,'defaulttextColor'),...
'FontSize',[],...
'Position',[0.525252525252525 0.743566176470588 0],...
'String','5 MPa',...
'ButtonDownFcn',mat{75},...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h29 = line(...
'Parent',h20,...
'Color',[0.455832079626395 0.789744115511128 0.645829808804391],...
'LineWidth',2,...
'XData',[0.0808080808080808 0.484848484848485],...
'YData',[0.743566176470588 0.743566176470588],...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'Tag','5 MPa');

h30 = line(...
'Parent',h20,...
'Color',[0.455832079626395 0.789744115511128 0.645829808804391],...
'LineStyle','none',...
'LineWidth',2,...
'XData',0.282828282828283,...
'YData',0.743566176470588,...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off');

appdata = [];
appdata.Listeners = mat{76};

h31 = text(...
'Parent',h20,...
'Units','normalized',...
'Color',get(0,'defaulttextColor'),...
'FontSize',[],...
'Position',[0.525252525252525 0.646139705882353 0],...
'String','10 MPa',...
'ButtonDownFcn',mat{77},...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h32 = line(...
'Parent',h20,...
'Color',[0.761635738765942 0.905831090606177 0.629864253393665],...
'LineWidth',2,...
'XData',[0.0808080808080808 0.484848484848485],...
'YData',[0.646139705882353 0.646139705882353],...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'Tag','10 MPa');

h33 = line(...
'Parent',h20,...
'Color',[0.761635738765942 0.905831090606177 0.629864253393665],...
'LineStyle','none',...
'LineWidth',2,...
'XData',0.282828282828283,...
'YData',0.646139705882353,...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off');

appdata = [];
appdata.Listeners = mat{78};

h34 = text(...
'Parent',h20,...
'Units','normalized',...
'Color',get(0,'defaulttextColor'),...
'FontSize',[],...
'Position',[0.525252525252525 0.548713235294118 0],...
'String','15 MPa',...
'ButtonDownFcn',mat{79},...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h35 = line(...
'Parent',h20,...
'Color',[0.92774517117992 0.958344124115124 0.644238952096614],...
'LineWidth',2,...
'XData',[0.0808080808080808 0.484848484848485],...
'YData',[0.548713235294118 0.548713235294118],...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'Tag','15 MPa');

h36 = line(...
'Parent',h20,...
'Color',[0.92774517117992 0.958344124115124 0.644238952096614],...
'LineStyle','none',...
'LineWidth',2,...
'XData',0.282828282828283,...
'YData',0.548713235294118,...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off');

appdata = [];
appdata.Listeners = mat{80};

h37 = text(...
'Parent',h20,...
'Units','normalized',...
'Color',get(0,'defaulttextColor'),...
'FontSize',[],...
'Position',[0.525252525252525 0.451286764705882 0],...
'String','20 MPa',...
'ButtonDownFcn',mat{81},...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h38 = line(...
'Parent',h20,...
'Color',[0.982025247227377 0.920603629584506 0.637231663608876],...
'LineWidth',2,...
'XData',[0.0808080808080808 0.484848484848485],...
'YData',[0.451286764705882 0.451286764705882],...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'Tag','20 MPa');

h39 = line(...
'Parent',h20,...
'Color',[0.982025247227377 0.920603629584506 0.637231663608876],...
'LineStyle','none',...
'LineWidth',2,...
'XData',0.282828282828283,...
'YData',0.451286764705882,...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off');

appdata = [];
appdata.Listeners = mat{82};

h40 = text(...
'Parent',h20,...
'Units','normalized',...
'Color',get(0,'defaulttextColor'),...
'FontSize',[],...
'Position',[0.525252525252525 0.353860294117647 0],...
'String','25 MPa',...
'ButtonDownFcn',mat{83},...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h41 = line(...
'Parent',h20,...
'Color',[0.994219317356572 0.758258112279141 0.431162294534286],...
'LineWidth',2,...
'XData',[0.0808080808080808 0.484848484848485],...
'YData',[0.353860294117647 0.353860294117647],...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'Tag','25 MPa');

h42 = line(...
'Parent',h20,...
'Color',[0.994219317356572 0.758258112279141 0.431162294534286],...
'LineStyle','none',...
'LineWidth',2,...
'XData',0.282828282828283,...
'YData',0.353860294117647,...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off');

appdata = [];
appdata.Listeners = mat{84};

h43 = text(...
'Parent',h20,...
'Units','normalized',...
'Color',get(0,'defaulttextColor'),...
'FontSize',[],...
'Position',[0.525252525252525 0.256433823529412 0],...
'String','30 MPa',...
'ButtonDownFcn',mat{85},...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h44 = line(...
'Parent',h20,...
'Color',[0.968399903171145 0.479865737670985 0.272320395922429],...
'LineWidth',2,...
'XData',[0.0808080808080808 0.484848484848485],...
'YData',[0.256433823529412 0.256433823529412],...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'Tag','30 MPa');

h45 = line(...
'Parent',h20,...
'Color',[0.968399903171145 0.479865737670985 0.272320395922429],...
'LineStyle','none',...
'LineWidth',2,...
'XData',0.282828282828283,...
'YData',0.256433823529412,...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off');

appdata = [];
appdata.Listeners = mat{86};

h46 = text(...
'Parent',h20,...
'Units','normalized',...
'Color',get(0,'defaulttextColor'),...
'FontSize',[],...
'Position',[0.525252525252525 0.159007352941177 0],...
'String','40 MPa',...
'ButtonDownFcn',mat{87},...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h47 = line(...
'Parent',h20,...
'Color',[0.852513916020359 0.2653890494876 0.308190107318648],...
'LineWidth',2,...
'XData',[0.0808080808080808 0.484848484848485],...
'YData',[0.159007352941177 0.159007352941177],...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'Tag','40 MPa');

h48 = line(...
'Parent',h20,...
'Color',[0.852513916020359 0.2653890494876 0.308190107318648],...
'LineStyle','none',...
'LineWidth',2,...
'XData',0.282828282828283,...
'YData',0.159007352941177,...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off');

appdata = [];
appdata.Listeners = mat{88};

h49 = text(...
'Parent',h20,...
'Units','normalized',...
'Color',get(0,'defaulttextColor'),...
'FontSize',[],...
'Position',[0.525252525252525 0.0615808823529413 0],...
'String','50 MPa',...
'ButtonDownFcn',mat{89},...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h50 = line(...
'Parent',h20,...
'Color',[0.619607843137255 0.00392156862745098 0.258823529411765],...
'LineWidth',2,...
'XData',[0.0808080808080808 0.484848484848485],...
'YData',[0.0615808823529413 0.0615808823529413],...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'Tag','50 MPa');

h51 = line(...
'Parent',h20,...
'Color',[0.619607843137255 0.00392156862745098 0.258823529411765],...
'LineStyle','none',...
'LineWidth',2,...
'XData',0.282828282828283,...
'YData',0.0615808823529413,...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off');

appdata = [];
appdata.CallbackObject = mat{90};

h52 = uicontextmenu(...
'Parent',h1,...
'HandleVisibility','off',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h53 = uimenu(...
'Parent',h52,...
'Callback',mat{91},...
'Label','Refresh',...
'HandleVisibility','off',...
'Tag','scribe:legend:refresh');

h54 = uimenu(...
'Parent',h52,...
'Callback',mat{92},...
'Label','Delete',...
'HandleVisibility','off',...
'Tag','scribe:legend:delete');

h55 = uimenu(...
'Parent',h52,...
'Callback',mat{93},...
'Label','Color ...',...
'Separator','on',...
'HandleVisibility','off',...
'Tag','scribe:legend:color');

h56 = uimenu(...
'Parent',h52,...
'Callback',mat{94},...
'Label','Edge Color ...',...
'HandleVisibility','off',...
'Tag','scribe:legend:edgecolor');

h57 = uimenu(...
'Parent',h52,...
'Callback',mat{95},...
'Label','Line Width',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth');

h58 = uimenu(...
'Parent',h57,...
'Callback',mat{96},...
'Label','0.5',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:0.5');

h59 = uimenu(...
'Parent',h57,...
'Callback',mat{97},...
'Label','1.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:1.0');

h60 = uimenu(...
'Parent',h57,...
'Callback',mat{98},...
'Label','2.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:2.0');

h61 = uimenu(...
'Parent',h57,...
'Callback',mat{99},...
'Label','3.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:3.0');

h62 = uimenu(...
'Parent',h57,...
'Callback',mat{100},...
'Label','4.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:4.0');

h63 = uimenu(...
'Parent',h57,...
'Callback',mat{101},...
'Label','5.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:5.0');

h64 = uimenu(...
'Parent',h57,...
'Callback',mat{102},...
'Label','6.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:6.0');

h65 = uimenu(...
'Parent',h57,...
'Callback',mat{103},...
'Label','7.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:7.0');

h66 = uimenu(...
'Parent',h57,...
'Callback',mat{104},...
'Label','8.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:8.0');

h67 = uimenu(...
'Parent',h57,...
'Callback',mat{105},...
'Label','9.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:9.0');

h68 = uimenu(...
'Parent',h57,...
'Callback',mat{106},...
'Label','10.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:10.0');

h69 = uimenu(...
'Parent',h57,...
'Callback',mat{107},...
'Label','11.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:11.0');

h70 = uimenu(...
'Parent',h57,...
'Callback',mat{108},...
'Label','12.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:12.0');

h71 = uimenu(...
'Parent',h52,...
'Callback',mat{109},...
'Label','Font ...',...
'HandleVisibility','off',...
'Tag','scribe:legend:font');

h72 = uimenu(...
'Parent',h52,...
'Callback',mat{110},...
'Label','Interpreter',...
'HandleVisibility','off',...
'Tag','scribe:legend:interpreter');

h73 = uimenu(...
'Parent',h72,...
'Callback',mat{111},...
'Label','latex',...
'HandleVisibility','off',...
'Tag','scribe:legend:interpreter:latex');

h74 = uimenu(...
'Parent',h72,...
'Callback',mat{112},...
'Label','tex',...
'HandleVisibility','off',...
'Tag','scribe:legend:interpreter:tex');

h75 = uimenu(...
'Parent',h72,...
'Callback',mat{113},...
'Label','none',...
'HandleVisibility','off',...
'Tag','scribe:legend:interpreter:none');

h76 = uimenu(...
'Parent',h52,...
'Callback',mat{114},...
'Label','Location',...
'HandleVisibility','off',...
'Tag','scribe:legend:location');

h77 = uimenu(...
'Parent',h76,...
'Callback',mat{115},...
'Label','Best',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:best');

h78 = uimenu(...
'Parent',h76,...
'Callback',mat{116},...
'Label','Inside North East',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:northeast');

h79 = uimenu(...
'Parent',h76,...
'Callback',mat{117},...
'Label','Outside North East',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:northeastoutside');

h80 = uimenu(...
'Parent',h76,...
'Callback',mat{118},...
'Label','Inside South East',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:southeast');

h81 = uimenu(...
'Parent',h76,...
'Callback',mat{119},...
'Label','Inside North West',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:northwest');

h82 = uimenu(...
'Parent',h76,...
'Callback',mat{120},...
'Label','Outside North West',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:northwestoutside');

h83 = uimenu(...
'Parent',h76,...
'Callback',mat{121},...
'Label','Inside South West',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:southwest');

h84 = uimenu(...
'Parent',h52,...
'Callback',mat{122},...
'Label','Orientation',...
'HandleVisibility','off',...
'Tag','scribe:legend:orientation');

h85 = uimenu(...
'Parent',h84,...
'Callback',mat{123},...
'Label','vertical',...
'HandleVisibility','off',...
'Tag','scribe:legend:orientation:vertical');

h86 = uimenu(...
'Parent',h84,...
'Callback',mat{124},...
'Label','horizontal',...
'HandleVisibility','off',...
'Tag','scribe:legend:orientation:horizontal');

h87 = uimenu(...
'Parent',h52,...
'Callback',mat{125},...
'Label','Show Property Editor',...
'Separator','on',...
'HandleVisibility','off',...
'Tag','scribe:legend:propedit');

h88 = uimenu(...
'Parent',h52,...
'Callback',mat{126},...
'Label','Show M-Code',...
'Separator','on',...
'HandleVisibility','off',...
'Tag','scribe:legend:mcode');

handles = [ h20 ];
set(handles, 'uicontextmenu', h52);



% --- Set application data first then calling the CreateFcn. 
function local_CreateFcn(hObject, eventdata, createfcn, appdata)

if ~isempty(appdata)
   names = fieldnames(appdata);
   for i=1:length(names)
       name = char(names(i));
       setappdata(hObject, name, getfield(appdata,name));
   end
end

if ~isempty(createfcn)
   if isa(createfcn,'function_handle')
       createfcn(hObject, eventdata);
   else
       eval(createfcn);
   end
end
