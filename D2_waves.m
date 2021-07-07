%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D-WAVES: CARTESIAN COORDINATES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
rng('shuffle');

% PARAMETERS %

% Differentials:
dt = 0.0005;
dx = sqrt(50/1000)/20;

% Nodes:
Nx = 50;
Ny = 50;  

x(:,1) = (0:Nx-1)*dx;
y(:,1) = (0:Ny-1)*dx; 

[X,Y] = meshgrid(x,y);

% Speed, dispersion, time:
E = 50;
p = 1000;

c = sqrt(E/1000);  
% c = sqrt(50/1000)*10^(6);  
r = (dt*c)/(dx);
D = -(E*dt/dx^(2))*10/p;     

d = D*dt;

frecuency = 5;
amplitude = 0.03;
T = 5000;
t(:,1)= (0:T-1)*dt;

% Position:
U = zeros(Nx,Ny,T);
U = reshape(U,Nx*Ny,T); 
Up = zeros(Nx);

% External force and sensor:
f_ext = zeros(Nx*Ny,1);
sensor_pos_kalman_x = round(Nx/2) + 10;
sensor_pos_kalman_y = round(Ny/2);

% MATRIX %

%Laplacian matrix:
A = laplace(Nx);

% UPDATE PREDICTION %
for k = 2:length(t)-1
    
    % External perturbation
%     f_ext(round((sqrt(Nx*Ny)/2 )+ Nx*Ny/2)) = amplitude*sin(frecuency*2*pi*t(k));
    f_ext(round((sqrt(Nx*Ny)/2 )+ Nx*Ny/2)) = amplitude*sin(frecuency*2*pi*t(k));
    
    % Bidimensional Waves Equation
	U(:,k+1) = ((r.^2).*A*(U(:,k)) - U(:,k)*((d)-2) - U(:,k-1) + (f_ext)*(dt^2))/(1-(d)); 
    % U(:,k+1) = ((r.^2).*A*(U(:,k)) - U(:,k)*(d-2) - U(:,k-1) + (f_ext)*(dt^2) + system_noise)/(1-d); 

    % Matrix remodelation 2D -> 3D
	Up(:,:,k+1) = reshape(U(:,k+1),Nx,Ny);
	sensor_kalman(:,k+1) = Up(sensor_pos_kalman_x ,sensor_pos_kalman_y,k);
        
end

% Noisy Kalman:
noise_per = 0.40;
random_numbers = randn(1,length(sensor_kalman))*(dt^2);
noisy_sensor = sensor_kalman + noise_per*random_numbers;

% per_noise = 0.1; % 10%
% base_noise = 0.1*sensor_kalman;
% random_numbers = randn(1,length(sensor_kalman));
% 
% for i = 1:T
%     base_noise(i) = base_noise(i)*random_numbers(i);
% end
% 
% noisy_sensor_1 = sensor_kalman+base_noise;

% PLOT %

% Kalman sensor:
figure(1);
plot1 = plot(t, sensor_kalman, '-', 'MarkerFaceColor',[0.49,0.18,0.56]) ;
% str = ' Time (t.u.) '; h=xlabel(str,'Interpreter','latex'); s=h.FontSize; h.FontSize=60;
% str = ' Position (p.u.) '; h=ylabel(str,'Interpreter','latex'); s=h.FontSize; h.FontSize=60;
str = ' $c_{f}/c_{r}$ '; h=xlabel(str,'Interpreter','latex'); s=h.FontSize; h.FontSize=60;
str = ' $d_{f}/d_{r}$ '; h=ylabel(str,'Interpreter','latex'); s=h.FontSize; h.FontSize=60;
ax = gca;
ax.FontSize = 25;
ylim([-4*10^(-6), 4*10^(-6)])
xlim([0, 2.5])
grid


% Noisy sensor:
figure(2);
plot1 = plot(t, noisy_sensor, 'MarkerFaceColor',[0.49,0.18,0.56]) ;
str = ' Time (t.u.) '; h=xlabel(str,'Interpreter','latex'); s=h.FontSize; h.FontSize=60;
str = ' Position (p.u.) '; h=ylabel(str,'Interpreter','latex'); s=h.FontSize; h.FontSize=60;
ax = gca;
ax.FontSize = 25;
ylim([-4*10^(-6), 4*10^(-6)])
xlim([0, 2.5])
grid

% Save data:
save(['vector_2d_12param.mat'], 'sensor_kalman')
save(['vector_2d_12param_noise.mat'], 'noisy_sensor')

% % ANIMATION PLOTS %
% 
% % Wave animation: 
% fh = figure(1);
% set(fh, 'Color', 'white'); 
% set(fh,'Position',[10 50 1350 610]);
% 
% video = VideoWriter('animation_fast.avi', 'Uncompressed AVI'); 
% video.FrameRate = 4;
% open(video)
%  for i = 10:1:length(t)-1  
%    Uplot(:,:) = Up(:,:,i);%Up(:,Nx/2  + 1,i);
% 
%       huhu = surf(x,y,Uplot','EdgeColor','none','LineStyle','none','FaceLighting','phong');
%      colormap jet;
% 
%       lighting flat
%       shading interp %%%%
%    view([90 90]);
%       colorbar('location','eastoutside');
%       colorbar
%      caxis([-max(sensor_kalman*15) max(sensor_kalman*15)]);
%            axis([min(x) max(x) min(y) max(y) -max(sensor_kalman*10) max(sensor_kalman*10)]);
% 
% %     xlabel('X axis');
%   
%    h=gca; 
% 
%  F = getframe(fh);
%     writeVideo(video,F)
%  end
% 
% close(video);