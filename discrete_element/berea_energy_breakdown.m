function h1 = berea_energy_breakdown()
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


load berea_energy_breakdown.mat


h1 = figure(...
'Color',[1 1 1],...
'Colormap',[0.346666666666667 0.536 0.690666666666667;0.915294117647059 0.28156862745098 0.287843137254902;0.44156862745098 0.749019607843137 0.432156862745098],...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'Position',get(0,'defaultfigurePosition'));

appdata = [];
appdata.SeriesBaseLine = mat{1};
appdata.barseriesXTick = mat{2};
appdata.PlotHoldStyle = mat{3};
appdata.legendInfoAffectedObject = mat{4};
appdata.LegendColorbarText = mat{5};
appdata.LegendColorbarInnerList = mat{6};
appdata.LegendColorbarOuterList = mat{7};
appdata.inLayout = mat{8};
appdata.LegendComputePosCache = mat{9};
appdata.LegendPeerHandle = mat{10};

h2 = axes(...
'Parent',h1,...
'Box','on',...
'CameraPosition',[45 3750 17.3205080756888],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'NextPlot','add',...
'XColor',get(0,'defaultaxesXColor'),...
'XLim',[-4.5 94.5],...
'XLimMode','manual',...
'XTick',mat{11},...
'XTickLabel',['0 ';'2 ';'5 ';'10';'15';'20';'25';'30';'40';'50'],...
'XTickLabelMode','manual',...
'XTickMode','manual',...
'YColor',get(0,'defaultaxesYColor'),...
'YLim',[0 7500],...
'YLimMode','manual',...
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
'FontSize',[],...
'FontWeight','normal',...
'HorizontalAlignment','center',...
'LineStyle','-',...
'LineWidth',0.5,...
'Margin',2,...
'Position',[44.8859447004608 7642.54385964912 1.00010919874411],...
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
'FontSize',[],...
'FontWeight','normal',...
'HorizontalAlignment','center',...
'LineStyle','-',...
'LineWidth',0.5,...
'Margin',2,...
'Position',[44.8859447004608 -515.350877192981 1.00010919874411],...
'Rotation',0,...
'String','Confining Pressure (MPa)',...
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
'FontSize',[],...
'FontWeight','normal',...
'HorizontalAlignment','center',...
'LineStyle','-',...
'LineWidth',0.5,...
'Margin',2,...
'Position',[-13.2822580645161 3717.1052631579 1.00010919874411],...
'Rotation',90,...
'String','Energy (Joules)',...
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
'FontSize',[],...
'FontWeight','normal',...
'HorizontalAlignment','right',...
'LineStyle','-',...
'LineWidth',0.5,...
'Margin',2,...
'Position',[-21.2661290322581 8168.85964912281 1.00010919874411],...
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
appdata.LegendLegendInfo = mat{12};
appdata.LegendLegendInfoStruct = mat{13};
appdata.LegendLegendType = mat{14};
appdata.legend_texthandle = mat{15};

h7 = specgraph.barseries(...
'Parent',h2,...
'DisplayName','Fracture Energy (W_{frac})',...
'XData',mat{16},...
'XDataMode','manual',...
'XDataSource',blanks(0),...
'YData',mat{17},...
'YDataSource',blanks(0),...
'OldSwitchProps',[],...
'OldSwitchVals',[],...
'Initialized',1,...
'Dirty','clean',...
'RefreshMode','auto',...
'BaseValue',0,...
'BaseLine',[],...
'BarLayout','stacked',...
'BarWidth',0.9,...
'Horizontal','off',...
'LineWidth',1,...
'EdgeColor',[0 0 0],...
'FaceColor','flat',...
'LineStyle','-',...
'ShowBaseLine','on',...
'CDataMapping','scaled',...
'BarPeers',[],...
'SwitchProps',{  'EdgeColor'; 'BarLayout'; 'BarWidth'; 'LineWidth'; 'LineStyle'; 'BaseValue'; 'ShowBaseLine'; 'DisplayName'; 'Visible'; 'FaceColor'; 'YDataSource'; 'XDataSource' },...
'BaseLineListener',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.legend_hgbehavior = mat{18};

h8 = specgraph.baseline(...
'Parent',h2,...
'Color',[0 0 0],...
'XData',[-4.5 94.5],...
'YData',[0 0],...
'HandleVisibility','off',...
'XLimInclude','off',...
'YLimInclude','off',...
'ZLimInclude','off',...
'BaseValue',0,...
'BaseValueMode','auto',...
'InternalSet',mat{19},...
'Orientation','X',...
'Listener',[],...
'AxesListener',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.LegendLegendInfo = mat{20};
appdata.LegendLegendInfoStruct = mat{21};
appdata.LegendLegendType = mat{22};
appdata.legend_texthandle = mat{23};

h9 = specgraph.barseries(...
'Parent',h2,...
'DisplayName','Internal Strain Energy (W_{int})',...
'XData',mat{24},...
'XDataMode','manual',...
'XDataSource',blanks(0),...
'YData',mat{25},...
'YDataSource',blanks(0),...
'OldSwitchProps',[],...
'OldSwitchVals',[],...
'Initialized',1,...
'Dirty','clean',...
'RefreshMode','auto',...
'BaseValue',0,...
'BaseLine',[],...
'BarLayout','stacked',...
'BarWidth',0.9,...
'Horizontal','off',...
'LineWidth',1,...
'EdgeColor',[0 0 0],...
'FaceColor','flat',...
'LineStyle','-',...
'ShowBaseLine','on',...
'CDataMapping','scaled',...
'BarPeers',[],...
'SwitchProps',{  'EdgeColor'; 'BarLayout'; 'BarWidth'; 'LineWidth'; 'LineStyle'; 'BaseValue'; 'ShowBaseLine'; 'DisplayName'; 'Visible'; 'FaceColor'; 'YDataSource'; 'XDataSource' },...
'BaseLineListener',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.LegendLegendInfo = mat{26};
appdata.LegendLegendInfoStruct = mat{27};
appdata.LegendLegendType = mat{28};
appdata.legend_texthandle = mat{29};

h10 = specgraph.barseries(...
'Parent',h2,...
'DisplayName','Inferred Frictional Energy (W_{fric})',...
'XData',mat{30},...
'XDataMode','manual',...
'XDataSource',blanks(0),...
'YData',mat{31},...
'YDataSource',blanks(0),...
'OldSwitchProps',[],...
'OldSwitchVals',[],...
'Initialized',1,...
'Dirty','clean',...
'RefreshMode','auto',...
'BaseValue',0,...
'BaseLine',[],...
'BarLayout','stacked',...
'BarWidth',0.9,...
'Horizontal','off',...
'LineWidth',1,...
'EdgeColor',[0 0 0],...
'FaceColor','flat',...
'LineStyle','-',...
'ShowBaseLine','on',...
'CDataMapping','scaled',...
'BarPeers',[],...
'SwitchProps',{  'EdgeColor'; 'BarLayout'; 'BarWidth'; 'LineWidth'; 'LineStyle'; 'BaseValue'; 'ShowBaseLine'; 'DisplayName'; 'Visible'; 'FaceColor'; 'YDataSource'; 'XDataSource' },...
'BaseLineListener',[],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h11 = text(...
'Parent',h2,...
'Color',get(0,'defaulttextColor'),...
'HandleVisibility','off',...
'Tag','LegendDeleteProxy',...
'Visible','off');

h12 = text(...
'Parent',h2,...
'Units','normalized',...
'Color',get(0,'defaulttextColor'),...
'HorizontalAlignment','center',...
'String',' ',...
'DeleteFcn','legendcolorbarlayout(get(gcbo,''Parent''),''remove'')',...
'HandleVisibility','off',...
'HitTest','off',...
'Tag','LegendColorbarLayout');

h13 = text(...
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
appdata.LegendTempText = mat{32};
appdata.MWBYPASS_title = mat{33};
appdata.MWBYPASS_xlabel = mat{34};
appdata.MWBYPASS_ylabel = mat{35};
appdata.MWBYPASS_zlabel = mat{36};
appdata.NonDataObject = mat{37};
appdata.PostDeserializeFcn = mat{38};
appdata.LegendOldSize = mat{39};

h14 = scribe.legend(...
'Parent',h1,...
'Position',[0.143154761904762 0.726587301587302 0.410714285714286 0.18015873015873],...
'Box','on',...
'CameraPosition',[0.5 0.5 17.3205080756888],...
'CameraPositionMode','auto',...
'CLim',[1 3],...
'CLimMode','manual',...
'Color',[1 1 1],...
'ColorOrder',mat{40},...
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
'ButtonDownFcn',mat{41},...
'Interruptible','off',...
'Tag','legend',...
'UserData',struct(...
    'PlotHandle', [], ...
    'legendpos', 2, ...
    'LegendPosition', [0.143154761904762 0.726587301587302 0.410714285714286 0.18015873015873], ...
    'LabelHandles', [], ...
    'handles', [], ...
    'lstrings', { {  'Fracture Energy (W_{frac})'; 'Internal Strain Energy (W_{int})'; 'Inferred Frictional Energy (W_{fric})' } }, ...
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
'String',{  'Fracture Energy (W_{frac})' 'Internal Strain Energy (W_{int})' 'Inferred Frictional Energy (W_{fric})' },...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h15 = text(...
'Parent',h14,...
'Units','points',...
'Color',get(0,'defaulttextColor'),...
'FontSize',[],...
'Margin',0.01,...
'String','Inferred Frictional Energy (W_{fric})',...
'HandleVisibility','off',...
'Tag','temphackytext',...
'Visible','off');

appdata = [];
appdata.Listeners = mat{42};

h16 = text(...
'Parent',h14,...
'Units','normalized',...
'Color',get(0,'defaulttextColor'),...
'FontSize',[],...
'Position',[0.226086956521739 0.812775330396476 0],...
'String','Fracture Energy (W_{frac})',...
'ButtonDownFcn',mat{43},...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h17 = hggroup(...
'Parent',h14,...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'Tag','Fracture Energy (W_{frac})');

h18 = patch(...
'Parent',h17,...
'FaceVertexCData',[1;1;1;1;1],...
'EdgeColor',get(0,'defaultpatchEdgeColor'),...
'FaceColor','flat',...
'Faces',[1 2 3 4 5],...
'LineWidth',1,...
'Vertices',[0.0347826086956522 0.695594713656388;0.0347826086956522 0.929955947136564;0.208695652173913 0.929955947136564;0.208695652173913 0.695594713656388;0.0347826086956522 0.695594713656388],...
'HitTest','off');

appdata = [];
appdata.Listeners = mat{44};

h19 = text(...
'Parent',h14,...
'Units','normalized',...
'Color',get(0,'defaulttextColor'),...
'FontSize',[],...
'Position',[0.226086956521739 0.5 0],...
'String','Internal Strain Energy (W_{int})',...
'ButtonDownFcn',mat{45},...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h20 = hggroup(...
'Parent',h14,...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'Tag','Internal Strain Energy (W_{int})');

h21 = patch(...
'Parent',h20,...
'FaceVertexCData',[2;2;2;2;2],...
'EdgeColor',get(0,'defaultpatchEdgeColor'),...
'FaceColor','flat',...
'Faces',[1 2 3 4 5],...
'LineWidth',1,...
'Vertices',[0.0347826086956522 0.382819383259912;0.0347826086956522 0.617180616740088;0.208695652173913 0.617180616740088;0.208695652173913 0.382819383259912;0.0347826086956522 0.382819383259912],...
'HitTest','off');

appdata = [];
appdata.Listeners = mat{46};

h22 = text(...
'Parent',h14,...
'Units','normalized',...
'Color',get(0,'defaulttextColor'),...
'FontSize',[],...
'Position',[0.226086956521739 0.187224669603524 0],...
'String','Inferred Frictional Energy (W_{fric})',...
'ButtonDownFcn',mat{47},...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h23 = hggroup(...
'Parent',h14,...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'Tag','Inferred Frictional Energy (W_{fric})');

h24 = patch(...
'Parent',h23,...
'FaceVertexCData',[3;3;3;3;3],...
'EdgeColor',get(0,'defaultpatchEdgeColor'),...
'FaceColor','flat',...
'Faces',[1 2 3 4 5],...
'LineWidth',1,...
'Vertices',[0.0347826086956522 0.0700440528634361;0.0347826086956522 0.304405286343612;0.208695652173913 0.304405286343612;0.208695652173913 0.0700440528634361;0.0347826086956522 0.0700440528634361],...
'HitTest','off');

appdata = [];
appdata.CallbackObject = mat{48};

h25 = uicontextmenu(...
'Parent',h1,...
'HandleVisibility','off',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h26 = uimenu(...
'Parent',h25,...
'Callback',mat{49},...
'Label','Refresh',...
'HandleVisibility','off',...
'Tag','scribe:legend:refresh');

h27 = uimenu(...
'Parent',h25,...
'Callback',mat{50},...
'Label','Delete',...
'HandleVisibility','off',...
'Tag','scribe:legend:delete');

h28 = uimenu(...
'Parent',h25,...
'Callback',mat{51},...
'Label','Color ...',...
'Separator','on',...
'HandleVisibility','off',...
'Tag','scribe:legend:color');

h29 = uimenu(...
'Parent',h25,...
'Callback',mat{52},...
'Label','Edge Color ...',...
'HandleVisibility','off',...
'Tag','scribe:legend:edgecolor');

h30 = uimenu(...
'Parent',h25,...
'Callback',mat{53},...
'Label','Line Width',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth');

h31 = uimenu(...
'Parent',h30,...
'Callback',mat{54},...
'Label','0.5',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:0.5');

h32 = uimenu(...
'Parent',h30,...
'Callback',mat{55},...
'Label','1.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:1.0');

h33 = uimenu(...
'Parent',h30,...
'Callback',mat{56},...
'Label','2.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:2.0');

h34 = uimenu(...
'Parent',h30,...
'Callback',mat{57},...
'Label','3.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:3.0');

h35 = uimenu(...
'Parent',h30,...
'Callback',mat{58},...
'Label','4.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:4.0');

h36 = uimenu(...
'Parent',h30,...
'Callback',mat{59},...
'Label','5.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:5.0');

h37 = uimenu(...
'Parent',h30,...
'Callback',mat{60},...
'Label','6.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:6.0');

h38 = uimenu(...
'Parent',h30,...
'Callback',mat{61},...
'Label','7.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:7.0');

h39 = uimenu(...
'Parent',h30,...
'Callback',mat{62},...
'Label','8.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:8.0');

h40 = uimenu(...
'Parent',h30,...
'Callback',mat{63},...
'Label','9.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:9.0');

h41 = uimenu(...
'Parent',h30,...
'Callback',mat{64},...
'Label','10.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:10.0');

h42 = uimenu(...
'Parent',h30,...
'Callback',mat{65},...
'Label','11.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:11.0');

h43 = uimenu(...
'Parent',h30,...
'Callback',mat{66},...
'Label','12.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:12.0');

h44 = uimenu(...
'Parent',h25,...
'Callback',mat{67},...
'Label','Font ...',...
'HandleVisibility','off',...
'Tag','scribe:legend:font');

h45 = uimenu(...
'Parent',h25,...
'Callback',mat{68},...
'Label','Interpreter',...
'HandleVisibility','off',...
'Tag','scribe:legend:interpreter');

h46 = uimenu(...
'Parent',h45,...
'Callback',mat{69},...
'Label','latex',...
'HandleVisibility','off',...
'Tag','scribe:legend:interpreter:latex');

h47 = uimenu(...
'Parent',h45,...
'Callback',mat{70},...
'Label','tex',...
'HandleVisibility','off',...
'Tag','scribe:legend:interpreter:tex');

h48 = uimenu(...
'Parent',h45,...
'Callback',mat{71},...
'Label','none',...
'HandleVisibility','off',...
'Tag','scribe:legend:interpreter:none');

h49 = uimenu(...
'Parent',h25,...
'Callback',mat{72},...
'Label','Location',...
'HandleVisibility','off',...
'Tag','scribe:legend:location');

h50 = uimenu(...
'Parent',h49,...
'Callback',mat{73},...
'Label','Best',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:best');

h51 = uimenu(...
'Parent',h49,...
'Callback',mat{74},...
'Label','Inside North East',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:northeast');

h52 = uimenu(...
'Parent',h49,...
'Callback',mat{75},...
'Label','Outside North East',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:northeastoutside');

h53 = uimenu(...
'Parent',h49,...
'Callback',mat{76},...
'Label','Inside South East',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:southeast');

h54 = uimenu(...
'Parent',h49,...
'Callback',mat{77},...
'Label','Inside North West',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:northwest');

h55 = uimenu(...
'Parent',h49,...
'Callback',mat{78},...
'Label','Outside North West',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:northwestoutside');

h56 = uimenu(...
'Parent',h49,...
'Callback',mat{79},...
'Label','Inside South West',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:southwest');

h57 = uimenu(...
'Parent',h25,...
'Callback',mat{80},...
'Label','Orientation',...
'HandleVisibility','off',...
'Tag','scribe:legend:orientation');

h58 = uimenu(...
'Parent',h57,...
'Callback',mat{81},...
'Label','vertical',...
'HandleVisibility','off',...
'Tag','scribe:legend:orientation:vertical');

h59 = uimenu(...
'Parent',h57,...
'Callback',mat{82},...
'Label','horizontal',...
'HandleVisibility','off',...
'Tag','scribe:legend:orientation:horizontal');

h60 = uimenu(...
'Parent',h25,...
'Callback',mat{83},...
'Label','Show Property Editor',...
'Separator','on',...
'HandleVisibility','off',...
'Tag','scribe:legend:propedit');

h61 = uimenu(...
'Parent',h25,...
'Callback',mat{84},...
'Label','Show M-Code',...
'Separator','on',...
'HandleVisibility','off',...
'Tag','scribe:legend:mcode');

handles = [ h14 ];
set(handles, 'uicontextmenu', h25);



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