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

% Last Modified by GUIDE v2.5 06-Feb-2016 17:37:05

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
%% 1.2 
% syms x1 x2 x3 u1 alf
% % �������������� ������ ����������� ��������
% F=[-x1+x2*x3+u1;
%     -x2-x1*x3-0.568*x3+0.1976;
%     5.4574*(x2-x3)-0.5434];
% % ���������� ����������� ����� 
% S = solve('-x1+x2*x3+u1=0','-x2-x1*x3-0.568*x3+0.1976=0','5.4574*(x2-x3)-0.5434=0','x1,x2,x3');
% Sx1=vpa(S.x1,4)
% Sx2=vpa(S.x2,4)
% Sx3=vpa(S.x3,4)
% 
% % ������������ � ����������� ����� ����������
% for i=1:numel(S.x1)
%     dF(i,1) = diff(F(i),x1);
%     dF(i,2) = diff(F(i),x2);
%     dF(i,3) = diff(F(i),x3);
% end
% % ������������������ �������
% A1 = vpa(subs(dF,{x1 x2 x3},{S.x1(1) S.x1(2) S.x1(3)}),4);
% A2 = vpa(subs(dF,{x1 x2 x3},{S.x2(1) S.x2(2) S.x2(3)}),4);
% A3 = vpa(subs(dF,{x1 x2 x3},{S.x3(1) S.x3(2) S.x3(3)}),4);
% % ������������������ �������� ������
% A1_p=vpa((det(alf*eye(numel(S.x1))-A1)),4)
% A2_p=vpa((det(alf*eye(numel(S.x1))-A2)),4)
% A3_p=vpa((det(alf*eye(numel(S.x1))-A3)),4)
% 
% 
% 
% 
% j=1;
% for u1i=-18:1:-3
%     L1(:,j)= solve(subs(A1_p,u1,u1i));
%     L2(:,j) = solve(subs(A2_p,u1,u1i));
%     L3 (:,j)= solve(subs(A3_p,u1,u1i)); 
% j=j+1;
% end
% % 
% 
% 
% axes(handles.axes1)
% plot(real(L1),imag(L1),'b*'),grid on
% xlabel({'�������� �������� ��� 1-��','������������������� ��������'})
% axes(handles.axes2)
% plot(real(L2),imag(L2),'b*'),grid on
% xlabel({'�������� �������� ��� 2-��','������������������� ��������'})
% axes(handles.axes3)
% plot(real(L3),imag(L3),'b*'),grid on
% xlabel({'�������� �������� ��� 3-��','������������������� ��������'})
% 
% 
% 
% u1=-15;
%  h = 0.01;
%  t=0:h:100;
% x_1(1)=0;
% x_2(1)=0;
% x_3(1)=0;
% for i=1:(length(t)-1)
%     x_1(i+1)=x_1(i)+h*(-x_1(i) + x_2(i)*x_3(i)+u1);
%     x_2(i+1)=x_2(i)+h*(-x_2(i)-x_1(i)*x_3(i)-.568*x_3(i)+0.1976);
%     x_3(i+1)=x_3(i)+h*(5.4574*(x_2(i)-x_3(i))-.5434);
% end
% axes(handles.axes4)
% plot3(x_1,x_2,x_3), grid on
% 
% 
% 
% T=50;
% n_=[-18:1:-3];
% x1=1;
% x2=1;
% x3=1;
% for j=1:length(n_)   
%     n=n_(j);
%     for k=1:T
%         x_1(1)=x1*(1+(10^-7)*rand); 
%         x_2(1)=x2*(1+(10^-7)*rand);
%         x_3(1)=x3*(1+(10^-7)*rand);
%         for i=1:(length(t)-1)           
%             x_1(i+1)=x_1(i)+h*(-x_1(i) + x_2(i)*x_3(i)+n);
%             x_2(i+1)=x_2(i)+h*(-x_2(i)-x_1(i)*x_3(i)-.568*x_3(i)+0.1976);
%             x_3(i+1)=x_3(i)+h*(5.4574*(x_2(i)-x_3(i))-.5434);
%         end
% %       yk(j,k)=x_1(i+1);
%           axes(handles.axes5)
% %         plot(n,yk(j,k),'k-*'), grid on, hold on, title('�������������� ���������')     
%      plot(n,x_1(i),'k-*'), grid on, hold on, title('�������������� ���������')    
%     end
% 
% 
% end
%       hold off



%% 1.6

syms x1 x2 x3 m
% �������� ���������� ��. ������� 2
F=[-10*x1+m*(x2-x3);
    12-x2+x1*x3;
    -x1*x2-x3+12];

F1_in=inline(F(1))
F2_in=inline(F(2))
F3_in=inline(F(3))

mu=100:10:250;

% ���������� ����������� ����� 
S = solve(F(1),F(2),F(3),'x1,x2,x3');
Sx1=vpa(S.x1,4)
Sx2=vpa(S.x2,4)
Sx3=vpa(S.x3,4)

J=jacobian(F,[x1 x2 x3]);

    J1 = subs(J,{x1 x2 x3},{S.x1(1), S.x2(1), S.x3(1)});
    J2 = subs(J,[x1 x2 x3],{S.x1(2), S.x2(2), S.x3(2)});
    J3 = subs(J,[x1 x2 x3],{S.x1(3), S.x2(3), S.x3(3)});
    

    j=1;
    for m_i=mu
    God1(:,j)=eig(double(subs(J1,m,m_i)));
    God2(:,j)=eig(double(subs(J2,m,m_i)));
    God3(:,j)=eig(double(subs(J3,m,m_i)));
    j=j+1;
    end

axes(handles.axes1)
plot(real(God1),imag(God1),'b*'),grid on
xlabel({'�������� �������� ��� 1-��','������������������� ��������'})
axes(handles.axes2)
plot(real(God2),imag(God2),'b*'),grid on
xlabel({'�������� �������� ��� 2-��','������������������� ��������'})
axes(handles.axes3)
plot(real(God3),imag(God3),'b*'),grid on
xlabel({'�������� �������� ��� 3-��','������������������� ��������'})


  mu_kr = min(mu(find(real(God3(3, :)) > 0)))  

  h=0.01;
  t=[0:h:100];
  
[T, X] = ode45(@odefun,t,[1;0;0], [], mu_kr);


axes(handles.axes4)
plot3(X(:, 1), X(:, 2), X(:, 3)); grid on



% NumberA=150;
% % ������� ��������� 
% A=linspace(100,250,NumberA);
%  
% N=1000;
% % ����� ����� �������, ������� ������ � ���������
% Num2Bif=10;
% b=10;
% c=12;
% for j=1:NumberA
%     a=A(j);
%     X(1,:)=0;
%     X(2,:)=12;
%     X(3,:)=1;
%     h=1e-3;
%     for i=1:N-1
%         X(1,i+1)=X(1,i)+h*(-b*X(1,i)+a*(X(2,i)-X(3,i)));
%         X(2,i+1)=X(2,i)+h*(c-X(2,i)+X(1,i)*X(3,i));
%         X(3,i+1)=X(3,i)+h*(-X(1,i)*X(2,i)-X(3,i)+c);
%     end
%     Bif(j,1:Num2Bif)=X(1,N-Num2Bif+1:N);
% end
%  
% axes(handles.axes5)
% plot(A,Bif,'k*-','LineWidth',2)
% grid on
% xlabel('\mu')
% ylabel('X_1')


function dxdt = odefun(t, x, mu)
 
dxdt = zeros(3, 1);
dxdt(1) = - 10 * x(1) + mu * (x(2) - x(3));
dxdt(2) = 12 - x(2) + x(1) * x(3);
dxdt(3) = - x(1) * x(2) - x(3) + 12;
