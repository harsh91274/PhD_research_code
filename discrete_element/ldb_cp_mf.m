function h1 = ldb_cp_mf()
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


load ldb_cp_mf.mat


h1 = figure(...
'Color',[1 1 1],...
'Colormap',[0.346666666666667 0.536 0.690666666666667;0.915294117647059 0.28156862745098 0.287843137254902],...
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
'CameraPosition',[45 1750 17.3205080756888],...
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
'YLim',[0 3500],...
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
'FontSize',10,...
'FontWeight','normal',...
'HorizontalAlignment','center',...
'LineStyle','-',...
'LineWidth',0.5,...
'Margin',2,...
'Position',[44.8859447004608 3566.52046783626 1.00010919874411],...
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
'Position',[44.8859447004608 -240.497076023392 1.00010919874411],...
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
'FontSize',10,...
'FontWeight','normal',...
'HorizontalAlignment','center',...
'LineStyle','-',...
'LineWidth',0.5,...
'Margin',2,...
'Position',[-13.2822580645161 1734.64912280702 1.00010919874411],...
'Rotation',90,...
'String','Number of Microfractures',...
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
'Position',[-21.2661290322581 3812.13450292398 1.00010919874411],...
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
'DisplayName','Shear',...
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
'DisplayName','Tensile',...
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

h10 = text(...
'Parent',h2,...
'Color',get(0,'defaulttextColor'),...
'HandleVisibility','off',...
'Tag','LegendDeleteProxy',...
'Visible','off');

h11 = text(...
'Parent',h2,...
'Units','normalized',...
'Color',get(0,'defaulttextColor'),...
'HorizontalAlignment','center',...
'String',' ',...
'DeleteFcn','legendcolorbarlayout(get(gcbo,''Parent''),''remove'')',...
'HandleVisibility','off',...
'HitTest','off',...
'Tag','LegendColorbarLayout');

h12 = text(...
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
appdata.LegendTempText = mat{26};
appdata.MWBYPASS_title = mat{27};
appdata.MWBYPASS_xlabel = mat{28};
appdata.MWBYPASS_ylabel = mat{29};
appdata.MWBYPASS_zlabel = mat{30};
appdata.NonDataObject = mat{31};
appdata.PostDeserializeFcn = mat{32};
appdata.LegendOldSize = mat{33};

h13 = scribe.legend(...
'Parent',h1,...
'Position',[0.143154761904762 0.811507936507937 0.176785714285714 0.0952380952380952],...
'Box','on',...
'CameraPosition',[0.5 0.5 17.3205080756888],...
'CameraPositionMode','auto',...
'CLim',[1 2],...
'CLimMode','manual',...
'Color',[1 1 1],...
'ColorOrder',mat{34},...
'DrawMode','fast',...
'FontSize',9,...
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
'ButtonDownFcn',mat{35},...
'Interruptible','off',...
'Tag','legend',...
'UserData',struct(...
    'PlotHandle', [], ...
    'legendpos', 2, ...
    'LegendPosition', [0.143154761904762 0.811507936507937 0.176785714285714 0.0952380952380952], ...
    'LabelHandles', [], ...
    'handles', [], ...
    'lstrings', { {  'Shear'; 'Tensile' } }, ...
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
'String',{  'Shear' 'Tensile' },...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h14 = text(...
'Parent',h13,...
'Units','points',...
'Color',get(0,'defaulttextColor'),...
'FontSize',9,...
'Margin',0.01,...
'String','Tensile',...
'HandleVisibility','off',...
'Tag','temphackytext',...
'Visible','off');

appdata = [];
appdata.Listeners = mat{36};

h15 = text(...
'Parent',h13,...
'Units','normalized',...
'Color',get(0,'defaulttextColor'),...
'FontSize',9,...
'Position',[0.525252525252525 0.720833333333333 0],...
'String','Shear',...
'ButtonDownFcn',mat{37},...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h16 = hggroup(...
'Parent',h13,...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'Tag','Shear');

h17 = patch(...
'Parent',h16,...
'FaceVertexCData',[1;1;1;1;1],...
'EdgeColor',get(0,'defaultpatchEdgeColor'),...
'FaceColor','flat',...
'Faces',[1 2 3 4 5],...
'LineWidth',1,...
'Vertices',[0.0808080808080808 0.559166666666667;0.0808080808080808 0.8825;0.484848484848485 0.8825;0.484848484848485 0.559166666666667;0.0808080808080808 0.559166666666667],...
'HitTest','off');

appdata = [];
appdata.Listeners = mat{38};

h18 = text(...
'Parent',h13,...
'Units','normalized',...
'Color',get(0,'defaulttextColor'),...
'FontSize',9,...
'Position',[0.525252525252525 0.279166666666667 0],...
'String','Tensile',...
'ButtonDownFcn',mat{39},...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h19 = hggroup(...
'Parent',h13,...
'HitTest','off',...
'Interruptible','off',...
'SelectionHighlight','off',...
'Tag','Tensile');

h20 = patch(...
'Parent',h19,...
'FaceVertexCData',[2;2;2;2;2],...
'EdgeColor',get(0,'defaultpatchEdgeColor'),...
'FaceColor','flat',...
'Faces',[1 2 3 4 5],...
'LineWidth',1,...
'Vertices',[0.0808080808080808 0.1175;0.0808080808080808 0.440833333333333;0.484848484848485 0.440833333333333;0.484848484848485 0.1175;0.0808080808080808 0.1175],...
'HitTest','off');

appdata = [];
appdata.CallbackObject = mat{40};

h21 = uicontextmenu(...
'Parent',h1,...
'HandleVisibility','off',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

h22 = uimenu(...
'Parent',h21,...
'Callback',mat{41},...
'Label','Refresh',...
'HandleVisibility','off',...
'Tag','scribe:legend:refresh');

h23 = uimenu(...
'Parent',h21,...
'Callback',mat{42},...
'Label','Delete',...
'HandleVisibility','off',...
'Tag','scribe:legend:delete');

h24 = uimenu(...
'Parent',h21,...
'Callback',mat{43},...
'Label','Color ...',...
'Separator','on',...
'HandleVisibility','off',...
'Tag','scribe:legend:color');

h25 = uimenu(...
'Parent',h21,...
'Callback',mat{44},...
'Label','Edge Color ...',...
'HandleVisibility','off',...
'Tag','scribe:legend:edgecolor');

h26 = uimenu(...
'Parent',h21,...
'Callback',mat{45},...
'Label','Line Width',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth');

h27 = uimenu(...
'Parent',h26,...
'Callback',mat{46},...
'Label','0.5',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:0.5');

h28 = uimenu(...
'Parent',h26,...
'Callback',mat{47},...
'Label','1.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:1.0');

h29 = uimenu(...
'Parent',h26,...
'Callback',mat{48},...
'Label','2.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:2.0');

h30 = uimenu(...
'Parent',h26,...
'Callback',mat{49},...
'Label','3.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:3.0');

h31 = uimenu(...
'Parent',h26,...
'Callback',mat{50},...
'Label','4.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:4.0');

h32 = uimenu(...
'Parent',h26,...
'Callback',mat{51},...
'Label','5.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:5.0');

h33 = uimenu(...
'Parent',h26,...
'Callback',mat{52},...
'Label','6.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:6.0');

h34 = uimenu(...
'Parent',h26,...
'Callback',mat{53},...
'Label','7.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:7.0');

h35 = uimenu(...
'Parent',h26,...
'Callback',mat{54},...
'Label','8.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:8.0');

h36 = uimenu(...
'Parent',h26,...
'Callback',mat{55},...
'Label','9.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:9.0');

h37 = uimenu(...
'Parent',h26,...
'Callback',mat{56},...
'Label','10.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:10.0');

h38 = uimenu(...
'Parent',h26,...
'Callback',mat{57},...
'Label','11.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:11.0');

h39 = uimenu(...
'Parent',h26,...
'Callback',mat{58},...
'Label','12.0',...
'HandleVisibility','off',...
'Tag','scribe:legend:linewidth:12.0');

h40 = uimenu(...
'Parent',h21,...
'Callback',mat{59},...
'Label','Font ...',...
'HandleVisibility','off',...
'Tag','scribe:legend:font');

h41 = uimenu(...
'Parent',h21,...
'Callback',mat{60},...
'Label','Interpreter',...
'HandleVisibility','off',...
'Tag','scribe:legend:interpreter');

h42 = uimenu(...
'Parent',h41,...
'Callback',mat{61},...
'Label','latex',...
'HandleVisibility','off',...
'Tag','scribe:legend:interpreter:latex');

h43 = uimenu(...
'Parent',h41,...
'Callback',mat{62},...
'Label','tex',...
'HandleVisibility','off',...
'Tag','scribe:legend:interpreter:tex');

h44 = uimenu(...
'Parent',h41,...
'Callback',mat{63},...
'Label','none',...
'HandleVisibility','off',...
'Tag','scribe:legend:interpreter:none');

h45 = uimenu(...
'Parent',h21,...
'Callback',mat{64},...
'Label','Location',...
'HandleVisibility','off',...
'Tag','scribe:legend:location');

h46 = uimenu(...
'Parent',h45,...
'Callback',mat{65},...
'Label','Best',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:best');

h47 = uimenu(...
'Parent',h45,...
'Callback',mat{66},...
'Label','Inside North East',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:northeast');

h48 = uimenu(...
'Parent',h45,...
'Callback',mat{67},...
'Label','Outside North East',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:northeastoutside');

h49 = uimenu(...
'Parent',h45,...
'Callback',mat{68},...
'Label','Inside South East',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:southeast');

h50 = uimenu(...
'Parent',h45,...
'Callback',mat{69},...
'Label','Inside North West',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:northwest');

h51 = uimenu(...
'Parent',h45,...
'Callback',mat{70},...
'Label','Outside North West',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:northwestoutside');

h52 = uimenu(...
'Parent',h45,...
'Callback',mat{71},...
'Label','Inside South West',...
'HandleVisibility','off',...
'Tag','scribe:legend:location:southwest');

h53 = uimenu(...
'Parent',h21,...
'Callback',mat{72},...
'Label','Orientation',...
'HandleVisibility','off',...
'Tag','scribe:legend:orientation');

h54 = uimenu(...
'Parent',h53,...
'Callback',mat{73},...
'Label','vertical',...
'HandleVisibility','off',...
'Tag','scribe:legend:orientation:vertical');

h55 = uimenu(...
'Parent',h53,...
'Callback',mat{74},...
'Label','horizontal',...
'HandleVisibility','off',...
'Tag','scribe:legend:orientation:horizontal');

h56 = uimenu(...
'Parent',h21,...
'Callback',mat{75},...
'Label','Show Property Editor',...
'Separator','on',...
'HandleVisibility','off',...
'Tag','scribe:legend:propedit');

h57 = uimenu(...
'Parent',h21,...
'Callback',mat{76},...
'Label','Show M-Code',...
'Separator','on',...
'HandleVisibility','off',...
'Tag','scribe:legend:mcode');

handles = [ h13 ];
set(handles, 'uicontextmenu', h21);



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