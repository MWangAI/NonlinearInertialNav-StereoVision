
close all
clc
global g e3 w LM Rc pL pR mI

N = 1; % numbers of landmark
LM = [0.8 0.6 0]';
% u = randn(3,1);
% u = u/norm(u);
% theta = 0.2*pi;
% Rc    = expm(theta*Skew(u));
Rc = eye(3);
pL = [0.1 0 0]';
pR = [0.1 0 0]';
mI = [0.8 0.6 0]';


%% Initialization
w = 1;
g  = 9.8;
e3 = [0 0 -1]';

% system
u = [0 1 0]';
u = u/norm(u);
theta = 0.0*pi;
R0    = expm(theta*Skew(u));
vR0   = reshape(R0,9,1);
p0    = 2*[1 0 0.5]';
v0    = 2*w*[0 1 0]';
x0     = [vR0;p0;v0];

% observer
u = [1 2 3]';
u = u/norm(u);
theta = 0.5*pi;
Rhat0 = expm(theta*Skew(u));
vRhat0= reshape(Rhat0,9,1);
vhat0 = 0*randn(3,1);
phat0 = 0*randn(3,1);
rhat0 = 0*randn(3,1);
xhat0 = [vRhat0;phat0;vhat0;rhat0];

temp  = 0.0*randn(9);
P     = 10*eye(length(xhat0)-9)+temp*temp';
vP    = P(:);%reshape(P,(length(P))^2,1);

X = [x0;xhat0;vP];

Tspan       =[0 20]; %
[Tout,Yout] = ode45(@odefun,Tspan,X); 

vRout    = Yout(:,1:9);
pout     = Yout(:,10:12);
vout     = Yout(:,13:15);
vRhatout = Yout(:,16:24);
phatout  = Yout(:,25:27);
vhatout  = Yout(:,28:30);
rhatout  = Yout(:,31:33);


error = zeros(4,length(Tout));
for i=1:length(Tout)
    vR    = vRout(i,:);
    p     = pout(i,:)';
    v     = vout(i,:)';
    vRhat = vRhatout(i,:);
    phat  = phatout(i,:)';
    vhat  = vhatout(i,:)';
    rhat  = rhatout(i,:)';
     
    R    = reshape(vR,3,3);
    Rhat = reshape(vRhat,3,3);
    
    error(i,1) = sqrt(trace(eye(3)-R*Rhat')/4);
    error(i,2) = norm(p- phat);
    error(i,3) = norm(v- vhat);
    error(i,4) = norm(g*e3- rhat);
end

figure
subplot(2,2,1)
plot(Tout,error(:,1),'linewidth',2)
grid on
xlabel('$t(s)$','interpreter','latex')
ylabel('$|\tilde{R}|_I$','interpreter','latex')
subplot(2,2,2)
plot(Tout,error(:,2),'linewidth',2)
grid on
xlabel('$t(s)$','interpreter','latex')
ylabel('$\|p-\hat{p}\|$','interpreter','latex')
subplot(2,2,3)
plot(Tout,error(:,3),'linewidth',2)
grid on
xlabel('$t(s)$','interpreter','latex')
ylabel('$\|v-\hat{v}\|$','interpreter','latex')
subplot(2,2,4)
plot(Tout,error(:,4),'linewidth',2)
grid on
xlabel('$t(s)$','interpreter','latex')
ylabel('$\|\mathsf{g}-\hat{\mathsf{g}}\|$','interpreter','latex')

% figure
% plot3(pout(:,1),pout(:,2),pout(:,3))
% plot3(phatout(:,1),phatout(:,2),phatout(:,3))
% zlim([0 3])


function DX=odefun(t,X)

global g e3 w LM Rc pL pR mI

vR    = X(1:9);
p     = X(10:12);
v     = X(13:15);
vRhat = X(16:24);
phat  = X(25:27);
vhat  = X(28:30);
rhat  = X(31:33);
vP    = X(34:end);


R    = reshape(vR,3,3);
Rhat = reshape(vRhat,3,3);
P    = reshape(vP,sqrt(length(vP)),sqrt(length(vP)));


%%  Output measurements 

% IMU measurements
Omega = 0.1*[0 0 1]';
a     = R'*(-2*w^2*[cos(w*t) sin(w*t) 0]'-g*e3);
mB    = R'*mI;

% bearing measurements
[~,col] = size(LM);
Xval = 0.00;
xL      = zeros(size(LM));
xR      = zeros(size(LM));
for i=1:col
    P_i      = LM(:,i);
    temp    = Rc'*(R'*(P_i-p)-pL);
    xL(:,i) = temp/norm(temp);        
    xL(:,i) = BearingNoise(xL(:,i),Xval); % Add noise
    temp    = Rc'*(R'*(P_i-p)-pR);
    xR(:,i) = temp/norm(temp);    
    xR(:,i) = BearingNoise(xR(:,i),Xval); % Add noise
end

% output
eL      = zeros(size(LM));
eR      = zeros(size(LM));
p1      = LM(:,1);
eL(:,1) = (eye(3)-xL(:,1)*xL(:,1)')*Rc'*(Rhat'*(p1-phat)-pL);
eR(:,1) = (eye(3)-xR(:,1)*xR(:,1)')*Rc'*(Rhat'*(p1-phat)-pR);  
y =  Rc*(eL(:,1)+eR(:,1));



%% LTV
A   = [-Skew(Omega) eye(3) zeros(3);
       zeros(3) -Skew(Omega) eye(3);
       zeros(3) zeros(3) -Skew(Omega)];
 
PI1 = Rc*(eye(3)-xL(:,1)*xL(:,1)'+eye(3)-xR(:,1)*xR(:,1)')*Rc';

C   = [PI1 zeros(3) zeros(3)];

[row,~] = size(C);
Q  = 1*eye(row);%inv(PI1);%
V  = 1000*eye(size(P));%diag([50*ones(1,3) 10*ones(1,3) 10*ones(1,3)]);%

DP  = (A*P+P*A'-P*C'*Q*C*P+V);
vDP = DP(:);

K  = P*C'*Q;

Kp = K(1:3,1:3); %10*eye(3);%
Kv = K(4:6,1:3); %1*eye(3);%
Kr = K(7:9,1:3); % 0.5*eye(3);%

%% Dynamics

% real system
DR  = R*Skew(Omega);
Dp  = v;
Dv  = g*e3 + R*a;
vDR = reshape(DR,9,1);
Dx  = [vDR;Dp;Dv];

% innovation terms     
rho=[1/g 1]';    
sigmaR = -rho(1)*cross(rhat,e3)-rho(2)*cross(Rhat*mB,mI);
sigmap = Skew(sigmaR)*(phat-p1)-Rhat*Kp*y;
sigmav = Skew(sigmaR)*vhat-Rhat*Kv*y;
sigmar = Skew(sigmaR)*rhat-Rhat*Kr*y;


% observer system
DRhat  = Rhat*Skew(Omega-Rhat'*sigmaR);
Dphat  = vhat-sigmap;
Dvhat  = rhat + Rhat*a-sigmav;
Drhat  = -sigmar;
vDRhat = reshape(DRhat,9,1);
Dxhat  = [vDRhat;Dphat;Dvhat;Drhat];

 

DX = [Dx;Dxhat;vDP];

end


