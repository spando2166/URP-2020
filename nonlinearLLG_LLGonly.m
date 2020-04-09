    
% This is a FDTD code solving LLG Equations of single spin unit
% Static demagnetization with excitation applied from quasistatic H-field
% Valid for magnetoquasistatic regime
% 3-d Cartesian coordinate system
%% Universal parameters
EP0 = 8.854e-12; % Epsilon, free space permittivity
MU0 = 4*pi*1.0E-7; % Mu, free space permeability
CC = 3.0*1.e8; % Speed of light in free space
eta0 = 120*pi;
gamma = 2*pi*2.8e6*10^4*MU0; % Unit is C/kg

% Define material properties
Ms_ga = 1750;
H0_oe = 100;

%units conversion
H0 = H0_oe/(4*pi*10^-3); % unit is A/m
Ms = Ms_ga/(4*pi*10^-3); % unit is A/m tesla, 4*pi*Ms = 1750 gauss -> 4*pi*Ms*10^-4 is tesla -> 4*pi*Ms*10^-4/MU0 is A/m
wm = 2*pi*2.8e6*Ms_ga;  % saturation magnetization
w0_0 = 2*pi*2.8e6*H0_oe;  % external biasing field


Nx=0;
Ny=1;
Nz=0;
wr = sqrt((w0_0+(Nx-Nz)*wm)*(w0_0+(Ny-Nz)*wm)); % FMR estimation based on Kittel's equation
fr=wr/2/pi;
delta_H = 10; % line width in the unit of Oesterd
alpha = delta_H*2.8e6*2*pi/w0_0/2;  % damping constant

%% Set up the time steps
DT  =  1/fr/(200*100); % define FDTD time-step size- careful of error accumulation as DT is increased
Tstep_end = round(70e-9/(DT));% total number of time steps = time duration/time resolution
Tstep_int = 1000; % number of time steps between adjacent record points reduce to save memory
Tstep_record = floor(Tstep_end/Tstep_int); % # of recorded data points
Time_axis = 0:Tstep_int:Tstep_end;

%% Time domain update

% Initialize field components;

%demagnetization terms
Nx=0;
Ny=1;
Nz=0;

Mx = Ms;
My = 0;
Mz = 0;

Hx = -Nx*Mx;
Hy = -Ny*My;
Hz = H0-Nz*Mz;

% Initialize field matrices
Hx_obs = zeros(1,Tstep_record);
Hy_obs = zeros(1,Tstep_record);
Hz_obs = zeros(1,Tstep_record);
Mx_obs = zeros(1,Tstep_record);
My_obs = zeros(1,Tstep_record);
Mz_obs = zeros(1,Tstep_record);

Prefactor1 = -gamma*DT;
Prefactor2 = Prefactor1*alpha/Ms;

for Tstepp = 1 : Tstep_end
    Tstep=Tstepp-1;
    TT = Tstep*DT;    
    
    %Applied H-field width demagnetization-needs update if applying
    %excitation
    %{
    Hx = -Nx*Mx;
    Hy= -Ny*My;
    Hz=H0-Nz*Mz;
    %}
    
    %LLG Equations
    deltaMx = Prefactor1 .* (My .* Hz - Mz .* Hy) + Prefactor2 .* (My .* (Mx .* Hy - My .* Hx)- Mz .* (Mz .* Hx - Mx .* Hz));
    deltaMy = Prefactor1 .* (Mz .* Hx - Mx .* Hz) + Prefactor2 .* (Mz .* (My .* Hz - Mz .* Hy)- Mx .* (Mx .* Hy - My .* Hx));
    deltaMz = Prefactor1 .* (Mx .* Hy - My .* Hx) + Prefactor2 .* (Mx .* (Mz .* Hx - Mx .* Hz)- My .* (My .* Hz - Mz .* Hy));
    
    %Update Magnetization
    Mx = Mx + deltaMx;
    My = My + deltaMy;
    Mz = Mz + deltaMz;
    
    %Renormalize magnetization to Ms
    mag = sqrt(Mx .* Mx + My .* My + Mz .* Mz);
    Mx = Mx ./ mag * Ms;
    My = My ./ mag * Ms;
    Mz = Mz ./ mag * Ms;
    
    if rem(Tstep, Tstep_int) == 0
        NN = floor(Tstep/Tstep_int);
        Mx_obs(NN+1) = Mx;
        My_obs(NN+1) = My;
        Mz_obs(NN+1) = Mz;
    end
end

figure
plot(Time_axis(1:length(Mx_obs))*DT, Mx_obs/(1e3/4/pi), 'k', 'Linewidth', 2)
hold on
plot(Time_axis(1:length(My_obs))*DT, My_obs/(1e3/4/pi), ':r', 'Linewidth', 4)

figure
% yyaxis right
plot(Time_axis(1:length(Mz_obs))*DT, Mz_obs/(1e3/4/pi), '-.b', 'Linewidth', 2)

legend ('Mx','My','Mz (bias direction)')
xlabel('Time [s]');
ylabel('4\piM [Gauss]');
% title( ['Magnetization magnitude, \alpha = ', num2str(alpha), ', X-direction RF Amplitude = Ho*' num2str(Hrf)]);

%% gif
h = figure;
filename = 'testAnimated.gif';
SpinCurve = animatedline('LineWidth',1);
% set(gca,'XLim',[-1000 1000],'YLim',[1450 1760],'ZLim',[-250 250])
xlabel('Mx, excitation direction');
ylabel('My');
zlabel('Mz, bias direction');
view(24,24);

for t3D = 1:5:length(Mx_obs)
    addpoints(SpinCurve,Mx_obs(t3D)/(1e3/4/pi),My_obs(t3D)/(1e3/4/pi),Mz_obs(t3D)/(1e3/4/pi));
    drawnow;
    %         pause(0.02)
    
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to the GIF File
    if t3D == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.01);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.01);
    end
    
end

%{
%% plot only
h = figure;
SpinCurve = animatedline('LineWidth',1);
% set(gca,'XLim',[-1000 1000],'YLim',[1450 1760],'ZLim',[-250 250])
view(40,24);

for t3D = 1:5:length(Mx_obs)
    addpoints(SpinCurve,Mx_obs(t3D)/(1e3/4/pi),My_obs(t3D)/(1e3/4/pi),Mz_obs(t3D)/(1e3/4/pi));
    drawnow;
end

%% avi

v=VideoWriter('SpinTip.avi');
v.FrameRate = 200;
%     v.Duration = 5;

open(v);
figure;
SpinCurve = animatedline('LineWidth',2);
% set(gca,'XLim',[-1000 1000],'YLim',[1450 1760],'ZLim',[-250 250])
view(13,24);

for t3D = 1:5:length(Mx_obs)
    addpoints(SpinCurve,Mx_obs(t3D)/(1e3/4/pi),My_obs(t3D)/(1e3/4/pi),Mz_obs(t3D)/(1e3/4/pi));
    drawnow;
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close (v)
%}
