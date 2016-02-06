function varargout = lab3_main(varargin)
% LAB3_MAIN MATLAB code for lab3_main.fig
%      LAB3_MAIN, by itself, creates a new LAB3_MAIN or raises the existing
%      singleton*.
%
%      H = LAB3_MAIN returns the handle to a new LAB3_MAIN or the handle to
%      the existing singleton*.
%
%      LAB3_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LAB3_MAIN.M with the given input arguments.
%
%      LAB3_MAIN('Property','Value',...) creates a new LAB3_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before lab3_main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to lab3_main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lab3_main

% Last Modified by GUIDE v2.5 06-Feb-2016 14:04:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @lab3_main_OpeningFcn, ...
                   'gui_OutputFcn',  @lab3_main_OutputFcn, ...
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


% --- Executes just before lab3_main is made visible.
function lab3_main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lab3_main (see VARARGIN)

% Choose default command line output for lab3_main
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes lab3_main wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = lab3_main_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%%  variant 33
% 1.2
% -14:-4 (1)
% u2=0.1976
% Mn=0.5434
% u1=-18 -3

% 1.6
% 109 157

% 1.11
% 0.08 0.4

% 2.2
% 0.44 0.77
%  2.6 
% 1.5
clc

%% 1.2 
syms x1 x2 x3 u1 alf
F=[-x1+x2*x3+u1;
    -x2-x1*x3-0.568*x3+0.1976;
    5.4574*(x2-x3)-0.5434];
S = solve('-x1+x2*x3+u1=0','-x2-x1*x3-0.568*x3+0.1976=0','5.4574*(x2-x3)-0.5434=0','x1,x2,x3')

for i=1:numel(S.x1)
    dF(i,1) = diff(F(i),x1);
    dF(i,2) = diff(F(i),x2);
    dF(i,3) = diff(F(i),x3);
end

A1 = vpa(subs(dF,{x1 x2 x3},{S.x1(1) S.x1(2) S.x1(3)}),4)
A2 = vpa(subs(dF,{x1 x2 x3},{S.x2(1) S.x2(2) S.x2(3)}),4)
A3 = vpa(subs(dF,{x1 x2 x3},{S.x3(1) S.x3(2) S.x3(3)}),4)

A1_p=vpa(det(alf*eye(numel(S.x1))-A1),4)
A2_p=vpa(det(alf*eye(numel(S.x1))-A2),4)
A3_p=vpa(det(alf*eye(numel(S.x1))-A3),4)
% 
for i=1:numel(S.x1)
for u1i=-18:-3

    L1 = solve(subs(A1_p,u1,u1i));
    L2 = solve(subs(A2_p,u1,u1i));
    L3 = solve(subs(A3_p,u1,u1i)); 
switch i
    case 1
        axes(handles.axes1)
        plot(real(L1(1)),imag(L1(1)),'*'),hold on
         plot(real(L1(2)),imag(L1(2)),'*'),hold on
          plot(real(L1(3)),imag(L1(3)),'*'),hold on
grid on
    case 2
        axes(handles.axes2)
        plot(real(L2(1)),imag(L2(1)),'*'),hold on
         plot(real(L2(2)),imag(L2(2)),'*'),hold on
          plot(real(L2(3)),imag(L2(3)),'*'),hold on
grid on
    case 3
        axes(handles.axes3)
        plot(real(L3(1)),imag(L3(1)),'*'),hold on
         plot(real(L3(2)),imag(L3(2)),'*'),hold on
          plot(real(L3(3)),imag(L3(3)),'*'),hold on
grid on
end
end
end