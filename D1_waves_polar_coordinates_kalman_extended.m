%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D-WAVES: POLAR COORDINATES %
% WITH KALMAN                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc

% 2D DATA %
y_sensor = struct2array(load('vector_2d_1param_noise.mat'));

% PARAMETERS %

% Differentials:
dr = sqrt(50/1000)/20;
dt = 0.0005;

% Nodes:
Nx = 50;

% Speed, dispersion, time:
E = 50;
p = 1000;

c_real = sqrt(E/1000);
c_estimated = 1.00 * c_real;

d_real = ((E)*(dt/dr^(2)))*10/p;
d_estimated = 1.25 * d_real;


r_real = (c_real*sqrt(dt))/dr;
r_estimated = (c_estimated*sqrt(dt))/dr;

frecuency = 5;
amplitude = 0.03; 

T = 5000;
t = (0:T-1)*dt;

% Position and velocity
u_real = zeros(1,Nx)';
v_real = zeros(1,Nx)';
X_real = [u_real; v_real];

u_estimated = zeros(1,Nx)';
v_estimated = zeros(1,Nx)';
X_estimated_EKF = [u_estimated; v_estimated; c_estimated; d_estimated];

% External force and sensor:
pos_force = Nx/2 ;
F_ext = zeros(2*Nx,1);

pos_sensor = pos_force + 10;
sensor_matrix = zeros(1,2*Nx);
sensor_matrix(pos_sensor) = 1;

sensor_matrix_EFK = zeros(1,2*Nx+2);
sensor_matrix_EFK(pos_sensor) = 1;

F_ext_eq = amplitude*sin(frecuency*2*pi*t); 

% MATRIX %

% Matrix A_h:
        A_h = ones(Nx);
        A_h = tril(A_h,+1); % returns the elements on and below the kth diagonal of A.
        A_h = triu(A_h,-1); % returns the elements on and above the kth diagonal of A.
        A_h = ((A_h - eye(Nx)));
        A_h = (-2)*eye(Nx) + A_h;
        
% Boundary conditions:
        A_h(2,1) = 2;
        A_h(Nx-1,Nx) = 2;
        
% Matrix A_b:
        A_b = eye(Nx);
        A_b = -1*ones(1,Nx-1);
        A_b = eye(Nx) + diag(A_b,1);
                
        for i = 1:(Nx)
            A_b(:,i) = A_b(:,i)/i;
        end
        
% Final matrix:
A_f = A_h + A_b;

% Function matrix:
A_real = [eye(Nx), eye(Nx)*(dt);
    (r_real^2)*A_f, (1-(d_real)*dt)*eye(Nx)]; 
          
A_estimada = [eye(Nx), eye(Nx)*(dt);
    (r_estimated^2)*A_f, (1-(d_estimated)*dt)*eye(Nx)];          
          
A_real_EKF = [A_real, zeros(2*Nx,2);
               zeros(2,2*Nx), eye(2)];         
          
A_estimada_EKF = [A_estimada, zeros(2*Nx,2);
                  zeros(2,2*Nx), eye(2)];

% Jacobiano:
M = [zeros(Nx,2);
     ((sqrt(c_estimated)*sqrt(dt))/dr)^(2)*2*A_f*u_estimated, -dt*eye(Nx)*v_estimated];  
 
I = eye(2);

O = zeros(2,2*Nx);

% F = [A_estimada, M;
%           O    ,      I    ];
 
% Matrix B: External force
B_force = zeros(2*Nx+2, 1); 
B_force(pos_force,1) = 1; 

% Matrix B: External noise -> we suppose 0
B_noise = zeros(2*Nx+2,2);

% Matrix B: discretised
B_dt = dt*[B_force, B_noise];
B_dt_ruido = B_dt(:,1:2); 
B_dt_force = B_dt(:,1);

% Simulation of random process noise
% noise_dv = diag([7 2]);
noise_dv = diag([7 2]);
noise_matrix = randn(length(t),2)*noise_dv;
noise_var = noise_dv^2;
input = [F_ext_eq' noise_matrix];

% Sensor
Csens = sensor_matrix_EFK;
n_sensor = size(sensor_matrix,1);

% We initialise variables to store graphs
grintcf = zeros(length(t),2);
post_estimate_var = 40^2*eye(2*Nx + 2); % Covariance/ Posteriori state estimate: initial estimate is very uncertain

noise_ratio_sensor = std(y_sensor);
noise_dv_sensor = noise_ratio_sensor;
noise_var_sensor = (noise_dv_sensor^2)*eye(n_sensor);

% UPDATE PREDICTION %
 for k = 1:length(t) - 1
     
    % External force
    F_ext(pos_force) = amplitude*sin(frecuency*2*pi*t(k));
    
    % Update matrix
    X_real(:,k+1) = A_real*X_real(:,k) + F_ext*(dt^2)*2*pi;

    % Update sensor data
    y(:,k+1) = sensor_matrix*X_real(:,k+1);
     
 end

% UPDATE KALMAN FILTER % 
for k=1:length(t)-1
    post_estimate_x_EKF = X_estimated_EKF;
    
    % PARAMETERS THAT ARE ESTIMATED
    u_estimated = X_estimated_EKF(1:Nx);
    v_estimated = X_estimated_EKF(Nx+1:2*Nx);
    c_estimated = X_estimated_EKF(end-1);
    d_estimated = X_estimated_EKF(end);
    
    r_estimated = (c_estimated*sqrt(dt))/dr;
    
    % MATRIX THAT ARE ESTIMATED
    A_estimada_EKF = [A_estimada, zeros(2*Nx,2);
                  zeros(2,2*Nx), eye(2)];
                                                   
    M = [zeros(Nx,2);
     ((sqrt(c_estimated)*sqrt(dt))/dr)^(2)*2*A_f*u_estimated, -dt*eye(Nx)*v_estimated];
         
    I = eye(2);
    
    O = zeros(2,2*Nx);
    
    F = [A_estimada, M;
        O    ,      I    ];
    
    % Predictor: variance of the a priori estimate
    priori_estimate_var = F*post_estimate_var*F' + B_dt_ruido*noise_var*B_dt_ruido';
    
    % Corrector: Kalman gain
    kalman_gain = (priori_estimate_var*Csens')*inv(Csens*priori_estimate_var*Csens'+noise_var_sensor);

    kalman_gain_pos_sensor(1,k+1) = kalman_gain(pos_sensor);
    kalman_gain_v_sensor(1,k+1) = kalman_gain(pos_sensor+Nx);
    kalman_gain_c(1,k+1) = kalman_gain(end-1);
    kalman_gain_d(1,k+1) = kalman_gain(end);
    
    % State estimate X: posteriori state estimate
%     post_estimate_x_EKF = X_estimated_EKF;
    X_estimated_EKF = X_estimated_EKF + kalman_gain*(y_sensor(k) - Csens*X_estimated_EKF);
    
    % A posteriori error covariance:
    post_estimate_var = priori_estimate_var - (kalman_gain*Csens*priori_estimate_var);
 
    % Hacemos simétrica si no lo es por redondeo:
    %esto es relevante en precisión finita -> hay otras implementaciones más rápidas computacionalmente
    post_estimate_var = (post_estimate_var+post_estimate_var')/2;
    
    % Predictor: Simulation of the next estimate state:
    X_estimated_EKF = (A_estimada_EKF*X_estimated_EKF + B_dt_force*F_ext_eq(k)*dt);

    % output_var = sensor_matrix_EKF*post_estimate_var*sensor_matrix_EKF';
    % grintcf(:,k+1) = 2.57583*sqrt(diag(output_var));
    % intervalo confianza -> 99% 
    % http://mathworld.wolfram.com/ConfidenceInterval.html

    c_matrix(1,k+1)= X_estimated_EKF(end-1);
    c_kk(1,k+1) = kalman_gain(end-1);
    d_matrix(1,k+1)= X_estimated_EKF(end);
    d_kk(1,k+1) = kalman_gain(end);
    ysr(:,k+1)= X_estimated_EKF;
end

% Error percentege
error_per = (c_real-c_matrix(end))/c_real;

% PLOT %
figure(1);
plot1 = plot(t, y, 'MarkerFaceColor',[0,0.51,0.255]) ;
str = ' Time (t.u.) '; h=xlabel(str,'Interpreter','latex'); s=h.FontSize; h.FontSize=60;
str = ' Position (p.u.) '; h=ylabel(str,'Interpreter','latex'); s=h.FontSize; h.FontSize=60;
ax = gca;
ax.FontSize = 25;
ylim([-4*10^(-6), 4*10^(-6)])
xlim([0, 2.5])
grid

figure(2);
plot1 = plot(t,y,'-',t,y_sensor,'-',t,ysr(pos_sensor,:),'--') ;
str = ' Time (t.u.) '; h=xlabel(str,'Interpreter','latex'); s=h.FontSize; h.FontSize=60;
str = ' Position (p.u.) '; h=ylabel(str,'Interpreter','latex'); s=h.FontSize; h.FontSize=60;
legend('Location','northwest')
h = legend('Model','Simulated data','Estimation with EKF','Interpreter','latex');
ax = gca;
ax.FontSize = 25;
ylim([-4*10^(-6), 4*10^(-6)])
xlim([0, 2.5])
grid


% PLOT 1D, 2D, Estimation with KFE
% figure(1);
% p = plot(t,y,'-g',t,y_sensor,'-r',t,ysr(pos_sensor,:),'b'); 
% set(p,'LineWidth',0.5)
% legend('Location','northwest')
% h = legend('Model','Simulated data','Estimation with EKF','Interpreter','latex');
% s = h.FontSize; 
% h.FontSize=9;
% h=xlabel('Time (t)','Interpreter','latex');
% s=h.FontSize;
% h.FontSize=15;
% h=ylabel('Position','Interpreter','latex');
% s=h.FontSize;
% h.FontSize=15;
% set(gcf,'position',[50,50,850,350])
% grid on
% title('Simulated data with Extended Kalman Filter');

% figure(2)
% plot(t, y,'^')
% str = '$$ \delta P\quad[Pa] $$';
% h=ylabel(str,'Interpreter','latex');
% s=h.FontSize;
% h.FontSize=18;
% str = '$$ time\quad[min] $$';
% h=xlabel(str,'Interpreter','latex');
% s=h.FontSize;
% h.FontSize=18;
% grid

% PLOT KALMAN GAIN
% plot(t, kalman_gain_pos_sensor,'r', t, kalman_gain_v_sensor,'b', t, kalman_gain_c, 'g', t, kalman_gain_d, 'm');legend('pos sensor','v sensor','c','d');
% xlabel('Time');
% ylabel('Kalman Gain');
% title('Kalman filter extended: Kalman Gain');

% PLOT KALMAN GAIN: POSITION AND VELOCITY  
% plot(t, kalman_gain_pos_sensor,'r', t, kalman_gain_v_sensor,'b');legend('pos_sensor','v_sensor');
% xlabel('Time');
% ylabel('Kalman Gain');
% title('Kalman filter extended: Kalman Gain');

% PLOT C AND D
figure(4)
c_vector = ones(1,T);
c_vector = c_real*c_vector;
p = plot(t(2:T), c_vector(2:T), 'r', t(2:T), c_matrix(2:T), 'b');
str = ' Time (t.u.) '; h=xlabel(str,'Interpreter','latex'); s=h.FontSize; h.FontSize=60;
str = ' Propagation Velocity $$c$$ (u.c.) '; h=ylabel(str,'Interpreter','latex'); s=h.FontSize; h.FontSize=60;
legend('Location','northwest')
h = legend('Estimated $$c$$','Real $$c$$','Interpreter','latex');
% ylim([min(c_matrix)*1.25 max(c_matrix)*1.25])
ylim([-0.0741  0.6571])
set(gcf,'position',[50,50,850,350])
ax = gca;
ax.FontSize = 25;
% set(gcf,'position',[900,400,450,350])
grid on
% % % 
figure(3)
d_vector = ones(1,T);
d_vector = d_real*d_vector;
p = plot(t(2:T), d_vector(2:T), 'r', t(2:T), d_matrix(2:T), 'b');
str = ' Time (t.u.) '; h=xlabel(str,'Interpreter','latex'); s=h.FontSize; h.FontSize=60;
str = ' Damping $$d$$ (u.d.) '; h=ylabel(str,'Interpreter','latex'); s=h.FontSize; h.FontSize=60;
legend('Location','northwest')
h = legend('Estimated d','Real d','Interpreter','latex');
ylim([-5.1532 16])
% ylim([min(d_matrix)*1.25 max(d_matrix)*1.25])
% set(gcf,'position',[50,50,850,350])
ax = gca;
ax.FontSize = 25;
grid on

c_final = c_matrix(end)/c_real; % c_real c_estimated
d_final = d_matrix(end)/d_real; % d_real d_estimated

% % PLOT DIFFERENCE BETWEEN REAL DATA AND ESTIMATION WITH KFE
% % plot(t,y-y_sensor,'g'); legend('Difference between real data and estimation with KFE');
% % xlabel('Time');
% % ylabel('Position');
% % title('1D+2D Waves: Sensor with Kalman filter extended');
% 
% % PLOT WITH KALMAN
% % plot(t,y,'-r',t,y_sensor,'g'); legend('1D','2D');
% % xlabel('Time');
% % ylabel('Position');
% % title('1D+2D Waves: Sensor with Kalman filter Extended');
% 
% plot(t, grintcf(1,:));