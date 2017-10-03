
function varargout = DragLineGUI(varargin)
% This demo is modeled after Doug Hull excellent demo @
% http://blogs.mathworks.com/pick/2008/05/27/advanced-matlab-capture-mouse-movement/
% but has been revised to work in a file created with GUIDE.
% The two lines can be dragged with the mouse and the positions can be passed to the command window by pressing the push button. 
% It also shows how to pass the GUI handles to the user defined callback functions in the same way GUIDE does.

% Edit the above text to modify the response to help DragLineGUI

% Last Modified by GUIDE v2.5 18-Jun-2013 12:35:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DragLineGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @DragLineGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before DragLineGUI is made visible.
function DragLineGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DragLineGUI (see VARARGIN)

cla(handles.axes1); % Clear axes

% This is how the handles are past to the callback function (the same way
% GUIDE defines the callback)
set(hObject, 'WindowButtonUpFcn', @(hObject,eventdata)DragLineGUI('stopDragFcn',hObject,eventdata,guidata(hObject)));
% http://www.mathworks.com/help/matlab/creating_guis/customizing-callbacks-in-guide.html
% Other ways I tried did not work:
%set(hObject, 'WindowButtonUpFcn', {'stopDragFcn',hObject,eventdata,guidata(hObject)});
%set(hObject, 'WindowButtonUpFcn', {@stopDragFcn, guidata(hObject)});
% http://www.mathworks.com/help/matlab/creating_plots/function-handle-callbacks.html

handles.MinFitInt = 0.5;
handles.MaxFitInt = 30;
axis(handles.axes1, [0.1 100 0 100]);
set(handles.axes1, 'XScale', 'log');
set(handles.axes1, 'box', 'on');

% Plot lines
handles.h_b_line = line([handles.MinFitInt handles.MinFitInt], [-2 102],...
    'LineWidth', 1, 'Color', 'b', 'ButtonDownFcn', @(hObject,eventdata)DragLineGUI('startDragB_Fcn',hObject,eventdata,guidata(hObject)));

handles.h_r_line = line([handles.MaxFitInt handles.MaxFitInt], [-2 102],...
    'LineWidth', 1, 'Color', 'r', 'ButtonDownFcn', @(hObject,eventdata)DragLineGUI('startDragR_Fcn',hObject,eventdata,guidata(hObject)));


% Choose default command line output for DragLineGUI
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
%assignin('base','handles',handles) % save handles to work space for debugging



% --- Outputs from this function are returned to the command line.
function varargout = DragLineGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
b = get(handles.h_b_line, 'XData');
sprintf('X coordinate of blue line: %5.3f\n',b(1))
r = get(handles.h_r_line, 'XData');
sprintf('X coordinate of red line: %5.3f\n',r(1))

function startDragB_Fcn(hObject, eventdata, handles)
        set(handles.figure1, 'WindowButtonMotionFcn', @(hObject,eventdata)DragLineGUI('draggingB_Fcn',hObject,eventdata,guidata(hObject))); 
    
function draggingB_Fcn(hObject, eventdata, handles)
        pt = get(handles.axes1, 'CurrentPoint');
        set(handles.h_b_line, 'XData', pt(1)*[1 1]);    
        
        
 
function startDragR_Fcn(hObject, eventdata, handles)
        set(handles.figure1, 'WindowButtonMotionFcn', @(hObject,eventdata)DragLineGUI('draggingR_Fcn',hObject,eventdata,guidata(hObject)));

function draggingR_Fcn(hObject, eventdata, handles)
        pt = get(handles.axes1, 'CurrentPoint');
        set(handles.h_r_line, 'XData', pt(1)*[1 1]);  
        
        
        
function stopDragFcn(hObject, eventdata, handles)
        set(handles.figure1, 'WindowButtonMotionFcn', '');
 
