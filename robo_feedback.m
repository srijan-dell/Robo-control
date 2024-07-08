% Feedback linearizing controller design for a two-link robot manipulator
% Aim is to design a feedback linearizing control law so that actual angular 
% position will track a desired (reference) angular position

t0=0; tf=10;  % Initial and final time
x0=[0.1 0 0 0]';  %  Initial condition
[t,x]=ode23(@roboct1,[0:0.01:10],x0);
plot(t,x(:,1),'r')
hold on
plot(t,x(:,3),'b')

xlabel('Time (sec)')
ylabel('Positions q1,q2')

function xdot = roboct1(t,x)
xdot=zeros(4,1);
% Compute desired trajectory to be tracked

period=2 ; amp1 =0.15 ; amp2 = 0.15 ;
fact=2*pi/period;
sinf=sin(fact*t) ;
cosf = cos(fact*t);
qd=[amp1*sinf  amp2*cosf]' ;
qdp=fact*[amp1*cosf  -amp2*sinf]' ;
qdpp=-fact^2*qd;

% PD feedback linearizing control law
% Link parameters
m1= 1.1; m2 = 1 ; a1 =1 ; a2=1.05; g = 9.81 ; 
% Controller gain selection
kp=100 ; kv=15;

% Define tracking errors

e= qd-[x(1) x(2)]'; % Desired - Actual position 
ep=qdp-[x(3) x(4)]'; % Desired - Actual velocity

% Compute M(q) and other terms of the dynamic equations

M11=(m1+m2)*a1^2+m2*a2^2+ 2*m2*a1*a2*cos(x(2));
M12=m2*a2^2+m2*a1*a2*cos(x(2));
M22=m2*a2^2;
N1=-m2*a1*a2*(2*x(3)*x(4) +x(4)^2)*sin(x(2));
N1=N1+(m1+m2)*g*a1*cos(x(1))+m2*g*a2*cos(x(1)+x(2));
N2=m2*a1*a2*x(3)^2*sin(x(2))+m2*g*a2*cos(x(1)+x(2));

% Control torques
s1=qdpp(1)+kv*ep(1)+kp*e(1);
s2=qdpp(2)+kv*ep(2)+kp*e(2);
tau1=M11*s1+M12*s2+N1;
tau2=M12*s1+M22*s2+N2;

% Inversion of mass matrix for deriving state equations
% M =[ M11 M12; M12 M22];

det =M11*M22-M12^2;
MI11=M22/det;
MI12=-M12/det;
MI22=M11/det;

% State equations of two link serial link manipulator
xdot(1)=x(3);
xdot(2)=x(4);
xdot(3)=MI11*(-N1+tau1)+MI12*(-N2+tau2);
xdot(4)=MI12*(-N1+tau1)+MI22*(N2+tau2);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


