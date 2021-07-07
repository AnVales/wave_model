%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D-WAVES: POLAR COORDINATES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc

% 2D DATA %
y_sensor = struct2array(load('vector_2d.mat'));

% PARAMETERS %

% Differentials:
dr = (sqrt(50/1000))/20;
dt = 0.0005;

% Nodes:
Nx = 50;

% Speed, dispersion, time:
c = sqrt(50/1000);
d = ((50)*(dt/dr^(2)));
r = (c*sqrt(dt))/dr;
T = 6000;
t = (0:T-1)*dt;
p = 1000;

amplitud = 0.03;
frecuency = 5;

% Position and velocity:
u = zeros(1,Nx)';
v = zeros(1,Nx)';
X = [u; v];

% External force and sensor:
pos_force = Nx/2 ;
F_ext = zeros(2*Nx,1);
pos_sensor = pos_force + 10;
sensor_matrix = zeros(1,2*Nx);
sensor_matrix(1, pos_sensor) = 1;

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
A = [eye(Nx), eye(Nx)*(dt);
              (r^2)*A_f, (1-(d/p)*dt)*eye(Nx)]; 

% UPDATE PREDICTION %
 for k = 1:length(t) - 1
     
    % External force
    F_ext(pos_force) = amplitud*sin(frecuency*2*pi*t(k));
    
    % Update matrix
    X(:,k+1) = A*X(:,k) + F_ext*(dt^2)*2*pi;

    % Update sensor data
    y(:,k+1) = sensor_matrix*X(:,k+1);
     
 end
 
% PLOTS %
 
% External force plot:
plot(t,X(pos_sensor,:));
xlabel('Time');
ylabel('Position');
title('1D Waves: External force');

% Comparation plot:
% plot(t,y,'-r',t,y_sensor,'-b'); legend('1D','2D with noise');
% xlabel('Time');
% ylabel('Position');
% title('1D+2D Waves: Sensor');