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

function dxdt = odefun12(t, x, u1)
 
dxdt = zeros(3, 1);
dxdt(1) =-x(1)+x(2)*x(3)+u1;
dxdt(2) =-x(2)-x(1)*x(3)+-0.568*x(3)+0.1976;
% dxdt(3) =5.4574*(x(2)-x(3))-0.5434*sign(x(3));
dxdt(3) =5.4574*(x(2)-x(3))+0.5434;


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
%% ========================        1.2  ===================================
clc

set(handles.pushbutton1,'String','����������� ...');
set(handles.pushbutton1,'Enable','inactive');
set(handles.uitable1,'RowName',{'������ �����','������ �����'});

syms x1 x2 x3 u1
% ��������� ����������� ���������
F=[-x1+x2*x3+u1;
    -x2-x1*x3+0.568*x3+0.1976;
    5.4574*(x2-x3)+0.5434];

F1_in=inline(F(1))
F2_in=inline(F(2))
F3_in=inline(F(3))

% �������������� ��������
u=-18:1:-3;

% ���������� ����������� ����� 
S = solve(F(1)==0,F(2)==0, F(3)==0,'x1,x2,x3');
Sx1=vpa(S.x1,4)
Sx2=vpa(S.x2,4)
Sx3=vpa(S.x3,4)

% ��������� ��������
J=jacobian(F,[x1 x2 x3])
% ����������� ����������� �����
    J1 = subs(J,[x1 x2 x3],{S.x1(1), S.x2(1), S.x3(1)});
    J2 = subs(J,[x1 x2 x3],{S.x1(2), S.x2(2), S.x3(2)});
    J3 = subs(J,[x1 x2 x3],{S.x1(3), S.x2(3), S.x3(3)});
% ���������� ����������� ��������
    j=1;
    for u1_i=u
    God1(:,j)=eig(double(subs(J1,u1,u1_i)));
    God2(:,j)=eig(double(subs(J2,u1,u1_i)));
    God3(:,j)=eig(double(subs(J3,u1,u1_i)));
    j=j+1;
    end
% ���������� �������� ����������
axes(handles.axes1)
plot(real(God1(3,:)),imag(God1(3,:)),'b*',real(God1(2,:)),imag(God1(2,:)),'r*',real(God1(1,:)),imag(God1(1,:)),'g*')
grid on
xlabel('Real'),ylabel('Image');
title({'�������� �������� ��� 1-��','������������������� ��������'})
axes(handles.axes2)
plot(real(God2(3,:)),imag(God2(3,:)),'b*',real(God2(2,:)),imag(God2(2,:)),'r*',real(God2(1,:)),imag(God2(1,:)),'g*')
grid on
xlabel('Real'),ylabel('Image');
title({'�������� �������� ��� 2-��','������������������� ��������'})
axes(handles.axes3)
plot(real(God3(3,:)),imag(God3(3,:)),'b*',real(God3(2,:)),imag(God3(2,:)),'r*',real(God3(1,:)),imag(God3(1,:)),'g*')
grid on
xlabel('Real'),ylabel('Image');
title({'�������� �������� ��� 3-��','������������������� ��������'})
% ������ ������������ �������� ��� ���������
u1_kr = min(u(find(real(God2(2,:)) < 0))) 
% % ========================���������� ���������� ��� ������������ ��������
% ��� ��������������
  h=0.01;
% �������� �������������
  T=0:h:50;
% ��������� ��������������
[Ti, Xi] = ode45(@odefun12,T,[1;0;0], [], u1_kr);
axes(handles.axes4)
plot3(Xi(:,1), Xi(:, 2), Xi(:, 3)); grid on
title(['A�������� ��� \mu=',num2str(u1_kr)])
% ====================================���������� �������������� ���������
Tt=20;
ub=-18:0.5:-3;
x1=Xi(fix(numel(Xi(:,1))/2),1);
x2=Xi(fix(numel(Xi(:,1))/2),2);
x3=Xi(fix(numel(Xi(:,1))/2),3);
for j=1:length(ub)    
tic
    for k=1:Tt
        x_1(1,k)=x1*(1+(10^-7)*rand); 
        x_2(1,k)=x2*(1+(10^-7)*rand);
        x_3(1,k)=x3*(1+(10^-7)*rand);

        for i=1:(length(T)-1)  
            x_1(i+1,k)=x_1(i,k)+h*(-x_1(i,k)+x_2(i,k)*x_3(i,k)+ub(j));
            x_2(i+1,k)=x_2(i,k)+h*(-x_2(i,k)-x_1(i,k)*x_3(i,k)+0.568*x3+0.1976);
%             x_3(i+1,k)=x_3(i,k)+h*(5.4574*(x_2(i,k)-x_3(i,k))-0.5434*sign(x_3(i,k)));
x_3(i+1,k)=x_3(i,k)+h*(5.4574*(x_2(i,k)-x_3(i,k))+0.5434);
        end

    end
       axes(handles.axes5)
        plot(ub(j),x_1(i+1,:),'k-*'), grid on, hold on  , title('�������������� ���������')  

end
hold off
%% ����������� ���������� �������� (1 ������) 
u1_ex=[-14 -4]
set(handles.uitable1,'ColumnName',{'-14','-4'})
    % ��� ��������������
  h=0.01;
% �������� �������������
  T=0:h:50;
% 
  bet=1.1;
for ii=1:numel(u1_ex)  
% ������� ����������
[Ti, Xi] = ode45(@odefun12,T,[1;0;0], [], u1_ex(ii));
% ����� �� ������� ���������       
T_=Ti(1:inv(h):end);
X_=Xi(1:inv(h):end,:); 

for k = 1 : T(end)
%         r0 = delta * sqrt(3);
r0=norm(X_(k, :)*bet-X_(k,:));
        xxx = ode45(@odefun12, [T_(k) : 0.01 : T_(k + 1)],X_(k, :)*bet,[], u1_ex(ii));
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

 for M=1:numel(u1_ex)
for i=1:N_int
    for j=1:N_var-1
        % �������� ������� (������� ����������)
       x_1=X(1,j);
       x_2=X(2,j);
       x_3=X(3,j);
       X(1,j+1)=x_1+h*(u1_ex(M)-x_1+x_2*x_3);
       X(2,j+1)=x_2+h*(-x_2+x_3*(7.1e1./1.25e2)-x_1*x_3+1.976e-1);
       X(3,j+1)=x_3+h*(x_2.*5.4574-x_3.*5.4574+5.434e-1);
 
        % ��������������� �������
    A=[-1 x_3 x_2;
        -x_3 -1 71/125-x_1;
        0 27287/5000 -27287/5000];
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
plot(real(God1(3,:)),imag(God1(3,:)),'b*',real(God1(2,:)),imag(God1(2,:)),'r*',real(God1(1,:)),imag(God1(1,:)),'g*')
grid on
xlabel('Real'),ylabel('Image');
title({'�������� �������� ��� 1-��','������������������� ��������'})
axes(handles.axes2)
plot(real(God2(3,:)),imag(God2(3,:)),'b*',real(God2(2,:)),imag(God2(2,:)),'r*',real(God2(1,:)),imag(God2(1,:)),'g*')
grid on
xlabel('Real'),ylabel('Image');
title({'�������� �������� ��� 2-��','������������������� ��������'})
axes(handles.axes3)
plot(real(God3(3,:)),imag(God3(3,:)),'b*',real(God3(2,:)),imag(God3(2,:)),'r*',real(God3(1,:)),imag(God3(1,:)),'g*')
grid on
xlabel('Real'),ylabel('Image');
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
hold off
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
