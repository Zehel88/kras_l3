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

% Last Modified by GUIDE v2.5 17-Feb-2016 17:15:45

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
% clc
format compact


function dxdt = odefun16(t, x, mu)
 
dxdt = zeros(3, 1);
dxdt(1) = - 10 * x(1) + mu * (x(2) - x(3));
dxdt(2) = 12 - x(2) + x(1) * x(3);
dxdt(3) = - x(1) * x(2) - x(3) + 12;


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
if get(handles.radiobutton1,'Value')==1
set(handles.radiobutton2,'Value',0);
set(handles.radiobutton3,'Value',0);
set(handles.radiobutton4,'Value',0);
set(handles.radiobutton5,'Value',0);
end

% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
if get(handles.radiobutton2,'Value')==1
set(handles.radiobutton1,'Value',0);
set(handles.radiobutton3,'Value',0);
set(handles.radiobutton4,'Value',0);
set(handles.radiobutton5,'Value',0);
end

% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
if get(handles.radiobutton3,'Value')==1
set(handles.radiobutton2,'Value',0);
set(handles.radiobutton1,'Value',0);
set(handles.radiobutton4,'Value',0);
set(handles.radiobutton5,'Value',0);
end

% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
if get(handles.radiobutton4,'Value')==1
set(handles.radiobutton2,'Value',0);
set(handles.radiobutton3,'Value',0);
set(handles.radiobutton1,'Value',0);
set(handles.radiobutton5,'Value',0);
end

% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
if get(handles.radiobutton5,'Value')==1
set(handles.radiobutton2,'Value',0);
set(handles.radiobutton3,'Value',0);
set(handles.radiobutton4,'Value',0);
set(handles.radiobutton1,'Value',0);
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% �������� ����� ��������� �������
active_var=find([get(handles.radiobutton1,'Value');
            get(handles.radiobutton2,'Value');
            get(handles.radiobutton3,'Value');
            get(handles.radiobutton4,'Value');
            get(handles.radiobutton5,'Value')] ==1);
        if isempty(active_var)==1
            errordlg('�� �� ������� ������� �������!');
        else
            
    switch active_var
        case 1
        case 2
%% ========================        1.6  ===================================
clc
set(handles.pushbutton1,'String','����������� ...');
set(handles.pushbutton1,'Enable','inactive');
set(handles.uitable1,'RowName',{'������ �����','������ �����'});

syms x1 x2 x3 m
% �������� ���������� ��. ������� 2
F=[-10*x1+m*(x2-x3);
    12-x2+x1*x3;
    -x1*x2-x3+12];

F1_in=inline(F(1))
F2_in=inline(F(2))
F3_in=inline(F(3))

% �������������� ��������
mu=100:5:250;

% ���������� ����������� ����� 
S = solve(F(1),F(2),F(3),'x1,x2,x3');
Sx1=vpa(S.x1,4)
Sx2=vpa(S.x2,4)
Sx3=vpa(S.x3,4)

% ��������� ��������
J=jacobian(F,[x1 x2 x3])
% ����������� ����������� �����
    J1 = subs(J,{x1 x2 x3},{S.x1(1), S.x2(1), S.x3(1)});
    J2 = subs(J,[x1 x2 x3],{S.x1(2), S.x2(2), S.x3(2)});
    J3 = subs(J,[x1 x2 x3],{S.x1(3), S.x2(3), S.x3(3)});
% ���������� ����������� ��������
    j=1;
    for m_i=mu
    God1(:,j)=eig(double(subs(J1,m,m_i)));
    God2(:,j)=eig(double(subs(J2,m,m_i)));
    God3(:,j)=eig(double(subs(J3,m,m_i)));
    j=j+1;
    end
% ���������� �������� ����������
axes(handles.axes1)
plot(real(God1),imag(God1),'b*'),grid on
title({'�������� �������� ��� 1-��','������������������� ��������'})
axes(handles.axes2)
plot(real(God2),imag(God2),'b*'),grid on
title({'�������� �������� ��� 2-��','������������������� ��������'})
axes(handles.axes3)
plot(real(God3),imag(God3),'b*'),grid on
title({'�������� �������� ��� 3-��','������������������� ��������'})
% ������ ������������ �������� ��� ���������
  mu_kr = min(mu(find(real(God3(3, :)) > 0)))  

% % ========================���������� ���������� ��� ������������ ��������
% ��� ��������������
  h=0.001;
% �������� �������������
  T=0:h:100;
% ��������� ��������������
[Ti, Xi] = ode45(@odefun16,T,[1;0;0], [], mu_kr);
axes(handles.axes4)
plot3(Xi(:,1), Xi(:, 2), Xi(:, 3)); grid on
title(['A�������� ��� \mu=',num2str(mu_kr)])
% ====================================���������� �������������� ���������
Tt=20;

x1=Xi(fix(numel(Xi(:,1))/2),1);
x2=Xi(fix(numel(Xi(:,1))/2),2);
x3=Xi(fix(numel(Xi(:,1))/2),3);
for j=1:length(mu)    
tic
    for k=1:Tt
        x_1(1,k)=x1*(1+(10^-7)*rand); 
        x_2(1,k)=x2*(1+(10^-7)*rand);
        x_3(1,k)=x3*(1+(10^-7)*rand);

        for i=1:(length(T)-1)      
            x_1(i+1,k)=x_1(i,k)+h*(-10*x_1(i,k)+(x_2(i,k)-x_3(i,k))*mu(j));
            x_2(i+1,k)=x_2(i,k)+h*(12-x_2(i,k)+x_1(i,k)*x_3(i,k));
            x_3(i+1,k)=x_3(i,k)+h*(-x_1(i,k)*x_2(i,k)-x_3(i,k)+12);
        end

    end
       axes(handles.axes5)
        plot(mu(j),x_1(i+1,:),'k-*'), grid on, hold on  , title('�������������� ���������')  

end

%% ����������� ���������� �������� (1 ������) 
mu_ex=[109 157]
set(handles.uitable1,'ColumnName',{'109','157'})
    % ��� ��������������
  h=0.001;
% �������� �������������
  T=0:h:50;
% 
  bet=1.1;
for ii=1:numel(mu_ex)  
% ������� ����������
[Ti, Xi] = ode45(@odefun16,T,[1;0;0], [], mu_ex(ii));
% ����� �� ������� ���������       
T_=Ti(1:inv(h):end);
X_=Xi(1:inv(h):end,:); 

for k = 1 : T(end)
%         r0 = delta * sqrt(3);
r0=norm(X_(k, :)*bet-X_(k,:));
        xxx = ode45(@odefun16, [T_(k) : 0.01 : T_(k + 1)],X_(k, :)*bet,[], mu_ex(ii));
        r1 = norm(xxx.y(:,end) - X_(k + 1, :)');
        lam(k)=log(r1/r0);
    end
 L1(ii)=1/T(end)*sum(lam); 
end
%%     ����������� ���������� �������� (2 ������)

tau=0.1;
N_var=tau/h;
N_int=T(end)/tau;
 
X(:,1)=[1;0;0];
 
% ��������� ����������
V(:,1)=rand(3,1);
V(:,1)=V(:,1)/norm(V(:,1));

 for M=1:numel(mu_ex)
for i=1:N_int
    for j=1:N_var-1
        % �������� ������� (������� ����������)
       x_1=X(1,j);
       x_2=X(2,j);
       x_3=X(3,j);
       X(1,j+1)=x_1+h*(-10*x_1+(x_2-x_3)*mu_ex(M));
       X(2,j+1)=x_2+h*(12-x_2+x_1*x_3);
       X(3,j+1)=x_3+h*(-x_1*x_2-x_3+12);
 
        % ��������������� �������
    A=[-10 mu_ex(M) -mu_ex(M);
        x_3 -1 x_1;
        -x_2 -x_1 -1];
        % ��������� � ���������
        V(:,j+1)=V(:,j)+h*A*V(:,j);
    end
    % ������������� ������ ��������� � �������� ���������� �������
    % ��� ���������� �������
    V(:,1)=V(:,j+1)/norm(V(:,j+1));
    X(:,1)=X(:,j+1);
    % ����� �������� 
    L(:,i)=norm(V(:,j+1));
    
end

L2(M)=(1/T(end))*sum(log(L));
 end
 
set(handles.uitable1,'ColumnWidth',{50,50});
set(handles.uitable1,'DaTa',[L1;L2]);

set(handles.pushbutton1,'String','���������');
set(handles.pushbutton1,'Enable','on');
        case 3
        case 4
        case 5
            
%     switch end
    end
% if empty
        end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
if get(handles.checkbox1,'Value')==1
    set(handles.checkbox1,'String','���');
else
    set(handles.checkbox1,'String','����');
end
