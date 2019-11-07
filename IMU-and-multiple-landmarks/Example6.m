
%ECC Conference paper
%Simulation with real trajectory and real IMU
%Virtual landmarks and landmark measurements
%Dataset from https://projects.asl.ethz.ch/datasets/doku.php?id=kmavvisualinertialdatasets


% close all
clear
clc


path        = 'datasetv2_01';
IMURead     = csvread([path '\imu0\data.csv'],1,0);
viconRead   = csvread([path '\vicon0\data.csv'],1,0);
Groundtruth = csvread([path '\state_groundtruth_estimate0\data.csv'],1,0);

timestart   = find(Groundtruth(1,1)==IMURead(:,1));
timeend     = find(Groundtruth(end,1)==IMURead(:,1));
IMURead     = IMURead(timestart:timeend,:);
bw          = mean(Groundtruth(:,12:14));
ba          = mean(Groundtruth(:,15:17));


% figure
% plot3(viconRead(:,2),viconRead(:,3),viconRead(:,4)), hold on
% plot3(Groundtruth(:,2),Groundtruth(:,3),Groundtruth(:,4))
% grid on

N = 5; % numbers of landmark
% Map = zeros(3,N);
% Map(1:2,:) = 4*randn(2,N);
Map = 10*[-0.5560   -0.1934    1.8198   -2.7144    0.8340;
        0.1306   1.4716   -0.2398   1.8803   -1.4323;
        0         0         0         0         0];

% u = randn(3,1);
% u = u/norm(u);
% theta = 0.2*pi;
% Rc    = expm(theta*Skew(u));
Rc = eye(3);
pL = [0 0 0]';
pR = [0.2 0 0]';
mI = [0.8 0.6 0]';


%% Initialization

g  = 9.8;
e3 = [0 0 1]';

p  = Groundtruth(1,2:4)';
R  = Q2R(Groundtruth(1,5:8)');
v  = Groundtruth(1,9:11)';
    
% observer
u = [1 0 0]';%randn(3,1);%
u = u/norm(u);
theta = 0.1*pi;
Rhat = expm(theta*Skew(u))*R; %R;%
phat = 0*p;
vhat = 0*v;
mhat = 0*Map(:,2:end);%randn(3,N-1);
rhat = -0*g*e3;


P     = 1*diag([1 1 1 1*ones(1,6+3*(N-1))]);




 

Tout  = IMURead(:,1)*1e-9;
Tout  = Tout - Tout(1,1); % set initial time as 0s
error = zeros(length(Tout));
L = zeros(size(Tout));


% figure
% subplot(1,2,1)
% plot(Tout,IMURead(:,2:4))
% xlim([0,100])
% xlabel('t(s)')
% ylabel('Gyroscope')
% subplot(1,2,2)
% plot(Tout,IMURead(:,5:7))
% xlim([0,100])
% xlabel('t(s)')
% ylabel('Accelerometer')
% set(gcf, 'Renderer', 'Painters');
% print('-depsc','IMURead.eps')



phatout = zeros(length(Tout),3); 
mhatout = zeros(length(Tout),3*(N-1));
p1      = Map(:,1); 
 
TTold = eye(6+3*N);
for k=1:length(Tout)
    if (k==1)
        dT = 0;
    else
        dT = (Tout(k,1)-Tout(k-1,1));  %nanoseconds to seconds
    end
    p  = Groundtruth(k,2:4)';
    R  = Q2R(Groundtruth(k,5:8)');
    v  = Groundtruth(k,9:11)';
 
    
    
    %%  Output measurements   
   
    % IMU measurements with bias corrected
    Omega = IMURead(k,2:4)' - 1*bw';
    a     = IMURead(k,5:7)' - 1*ba';
    
    
    % bearing measurements
    [~,N] = size(Map);
    Xval = 0.005;  % bearing measurement noise
    xL   = zeros(size(Map));
    xR   = zeros(size(Map));
    eL   = zeros(size(Map));
    eR   = zeros(size(Map));
    PI   = zeros(3,3*N);
    for i=1:N
        P_i      = Map(:,i);
        temp    = Rc'*(R'*(P_i-p)-pL);
        xL(:,i) = temp/norm(temp);        
        xL(:,i) = BearingNoise(xL(:,i),Xval); % Add noise
        temp    = Rc'*(R'*(P_i-p)-pR);
        xR(:,i) = temp/norm(temp);    
        xR(:,i) = BearingNoise(xR(:,i),Xval); % Add noise
    end 
    
    % virtual output
    p1      = Map(:,1);
   
    
    eL(:,1) = (eye(3)-xL(:,1)*xL(:,1)')*Rc'*(Rhat'*(p1-phat)-pL);
    eR(:,1) = (eye(3)-xR(:,1)*xR(:,1)')*Rc'*(Rhat'*(p1-phat)-pR);  
    y(1:3,1) =  Rc*(eL(:,1)+eR(:,1));
    PI(:,1:3) = Rc*(eye(3)-xL(:,1)*xL(:,1)'+eye(3)-xR(:,1)*xR(:,1)')*Rc';
    
%      p1 = mean(Map(:,2:end),2);
    
    for i=2:N      
        i3        = 3*(i-1);
        
        phati   = mhat(:,i-1);
        eL(:,i) = (eye(3)-xL(:,i)*xL(:,i)')*Rc'*(Rhat'*(phati-phat)-pL);
        eR(:,i) = (eye(3)-xR(:,i)*xR(:,i)')*Rc'*(Rhat'*(phati-phat)-pR); 
        temp    = Rc*((eye(3)-xL(:,i)*xL(:,i)')+(eye(3)-xR(:,i)*xR(:,i)'))*Rc';
        PI(:,i3+1:i3+3) = temp;
        y(i3+1:i3+3,1) = Rc*(eL(:,i)+eR(:,i));

 
    end

    %% LTV
    TempA = [-Skew(Omega) eye(3) zeros(3);
           zeros(3) -Skew(Omega) eye(3);
           zeros(3) zeros(3) -Skew(Omega)];
    Temp  = kron(eye(N-1),-Skew(Omega));
    A = [TempA zeros(9,3*(N-1));zeros(3*(N-1),9) Temp];
     
    Temp  = mat2cell(PI(:,4:end), 3, 3*ones(N-1,1));
    TempPI = blkdiag(Temp{:}); 
    
    C   = [PI(:,1:3) zeros(3) zeros(3) zeros(3,3*(N-1)); 
           PI(:,4:end)' zeros(3*(N-1),6) -TempPI]; 

    [row,~] = size(C); 
    
%     V  = 10*diag([1*ones(1,6) 0.1*ones(1,3+3*(N-1))]);%   5*eye(size(P));% 
     Q = 1000*eye(row);
     V  = 10*diag([2*ones(1,3) 2*ones(1,3) 1*ones(1,3) 0.1*ones(1,3*(N-1))]);%   5*eye(size(P)
    
% %     Approach 1 works

%     DP  = (A*P+P*A'-P*C'*Q*C*P+V);
%     P = P + DP*dT;
%     K = P*C'*Q;
%     
%     Approach 2 works (Runge-Kutta methods 4)
% Q = 0.001*eye(row);
% V = 0.001*eye(6+3*N);
%     DP1 = (A*P+P*A'-P*C'*Q*C*P+V);
%     P1  = P + DP1*dT/2;
%     DP2 = (A*P1+P1*A'-P1*C'*Q*C*P1+V);
%     P2  = P + DP2*dT/2;
%     DP3 = (A*P2+P2*A'-P2*C'*Q*C*P2+V);
%     P3  = P + DP3*dT;
%     DP4 = (A*P3+P3*A'-P3*C'*Q*C*P3+V);
%     P   = P + dT*(DP1+2*DP2+2*DP3+DP4)/6;
%     K   = P*C'*Q;


%     Approach 3 works
%      Ad = expm(A*dT); %discretized A matrix     
%      P  = Ad*P*Ad'+V;   
%      K  = P*C'/(C*P*C'+Q);
%      P  = (eye(size(P))-K*C)*P;
%      eig(P)

%       Approach 4 works
%      Abar = zeros(size(A));
%      Abar(1:3,4:6) = eye(3);
%      Abar(4:6,7:9) = eye(3);
%      TT = blkdiag(R,R,R,R,R,R,R);
%      Ad = TT'*expm(Abar*dT)*TTold;
%      TTold = TT;

% dicretized CRE
     Ad = expm(A*dT);
     P  = Ad*P*Ad'+V;
     K  = P*C'/(C*P*C'+eye(row)/Q);
     P  = (eye(size(P))-K*C)*P;

    

    Kp = 1*K(1:3,:); %10*eye(3);%
    Kv = 1*K(4:6,:); %1*eye(3);%
    Kr = 1*K(7:9,:); % 0.5*eye(3);%
    KI = 1*K(10:end,:);

    %% Dynamics 
   
    % innovation terms 
    sigmaR = zeros(3,1);
    for i=2:N
         sigmaR = sigmaR -(1/(N-1))*cross((mhat(:,i-1)-p1),(Map(:,i)-p1));
    end
     Kn = eye(N-1)/(N-1);
     M  = (Map(:,2:end)-p1)*Kn*(Map(:,2:end)-p1)';
     Mbar  = (trace(M)*eye(3)-M)/2;
     kR    = 0.4/max(eig(Mbar));
     sigmaR = kR*sigmaR;
    
    
%     sigmaR = - cross(Rhat*R'*g*e3,e3)/g-cross(Rhat*R'*mI,mI);
    
    sigmap = Skew(sigmaR)*(phat-p1)-Rhat*Kp*y;
    sigmav = Skew(sigmaR)*vhat-Rhat*Kv*y;
    sigmar = Skew(sigmaR)*rhat-Rhat*Kr*y;
    
    
    
    % observer system
    Omegahat = Omega*dT-Rhat'*sigmaR;
    Rhat  = Rhat*expm(Skew(Omegahat));
    phat  = phat + (vhat*dT-sigmap);
    vhat  = vhat + (rhat*dT + Rhat*a*dT-sigmav);
    rhat  = rhat + (-sigmar);
    for i=2:N
        i3 = 3*(i-2);
        Ki = KI(i3+1:i3+3,:);
        sigmapi = Skew(sigmaR)*(mhat(:,i-1)-p1) - Rhat*Ki*y;
        mhat(:,i-1) = mhat(:,i-1) + (-sigmapi); %Map(:,i);%
    end
     
    phatout(k,:) = phat';
%     vhatout(k,:) = vhat';
%     rhatout(k,:) = rhat'; 
    mhatout(k,:)=reshape(mhat,1,3*(N-1));
    %% estimation errors    
    error(k,1) = sqrt(trace(eye(3)-R*Rhat')/4);% max(eig(P));%
    error(k,2) = norm(p-phat);
    error(k,3) = norm(v-vhat);
    error(k,4) = norm(-g*e3-rhat); 
    
end

% eig(P)
clearvars -except Tout error  L viconRead Groundtruth phatout mhatout Map  

figure
hold on
index = 1:2:length(phatout);
% plot3(viconRead(:,2),viconRead(:,3),viconRead(:,4)), hold on
plot3(Groundtruth(index,2),Groundtruth(index,3),Groundtruth(index,4),'linewidth',1.0)
plot3(phatout(index,1),phatout(index,2),phatout(index,3),'linewidth',1.0)
for i=1:size(Map,2)-1
    i3= 3*(i-1);
    plot3(mhatout(index,i3+1),mhatout(index,i3+2),mhatout(index,i3+3),'c-.','linewidth',1)
    plot3(Map(1,i+1),Map(2,i+1),Map(3,i+1),'k*','linewidth',1.0)
end
legend('Ground truth' ,'Estimated state','Estimated landmarks','Landmarks')
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
grid on
 
% set(gcf, 'Renderer', 'Painters');
% print('-depsc','example6p1.eps')


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
% subplot(2,2,1)
% plot(Tout,error(:,1),'linewidth',1.5)
% grid on
% xlabel('$t(s)$','interpreter','latex')
% ylabel('$|\tilde{R}|_I$','interpreter','latex')
% subplot(2,2,2)
% plot(Tout,error(:,2),'linewidth',1.5)
% grid on
% xlabel('$t(s)$','interpreter','latex')
% ylabel('$\|\tilde{p}\|$','interpreter','latex')
% subplot(2,2,3)
% plot(Tout,error(:,3),'linewidth',1.5)
% grid on
% xlabel('$t(s)$','interpreter','latex')
% ylabel('$\|\tilde{v}\|$','interpreter','latex')
% subplot(2,2,4)
% plot(Tout,error(:,4),'linewidth',1.5)
% grid on
% xlabel('$t(s)$','interpreter','latex')
% ylabel('$\|\tilde{r}\|$','interpreter','latex')
% 
% figure
% % drawnow
% plot(Tout, L,'linewidth',2), hold on
% grid on
% xlabel('$t(s)$','interpreter','latex')
% ylabel('$\mathcal{L}$','interpreter','latex')

% figure
% subplot(3,1,1)
% plot(Tout,Groundtruth(:,2),Tout,phatout(:,1),'linewidth',1.5);
% grid on
% ylabel('x(m)')
% subplot(3,1,2)
% plot(Tout,Groundtruth(:,3),Tout,phatout(:,2),'linewidth',1.5);
% grid on
% ylabel('y(m)')
% subplot(3,1,3)
% plot(Tout,Groundtruth(:,4),Tout,phatout(:,3),'linewidth',1.5);
% grid on
% ylabel('z(m)')
% legend('Ground truth' ,'Estimated')
% xlabel('$t(s)$','interpreter','latex')
% set(gcf, 'Renderer', 'Painters');
% print('-depsc','example6p2.eps')
% 
% figure
% subplot(3,1,1)
% plot(Tout,Groundtruth(:,9),Tout,vhatout(:,1),'linewidth',1.5);
% grid on
% ylabel('x(m/s)')
% subplot(3,1,2)
% plot(Tout,Groundtruth(:,10),Tout,vhatout(:,2),'linewidth',1.5);
% grid on
% legend('Ground truth' ,'Estimated')
% 
% ylabel('y(m/s)')
% subplot(3,1,3)
% plot(Tout,Groundtruth(:,11),Tout,vhatout(:,3),'linewidth',1.5);
% grid on
% ylabel('z(m/s)')
% 
% xlabel('$t(s)$','interpreter','latex')
% 
% set(gcf, 'Renderer', 'Painters');
% print('-depsc','example6p3.eps')
% 
% figure
% subplot(2,1,1)
% plot(Tout,rhatout(:,1),'linewidth',1.5), hold on
% plot(Tout,rhatout(:,2),'linewidth',1.5);
% plot(Tout,rhatout(:,3),'linewidth',1.5);
% grid on
% legend('$\hat{r}_1$' ,'$\hat{r}_2$', '$\hat{r}_3$','interpreter','latex')
% xlabel('$t(s)$','interpreter','latex')
% ylabel('$\hat{r}$','interpreter','latex')
% 
% subplot(2,1,2)
% plot(Tout,error(:,1),'linewidth',1.5)
% grid on
% xlabel('$t(s)$','interpreter','latex')
% ylabel('$|\tilde{R}|_I$','interpreter','latex')

  



