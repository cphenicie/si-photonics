% MATLAB script to plot the MZI transfer function.
% by Lukas Chrostowski
% user must configure several variables below.
        
% specify the wavelength range of interest.
lambda_min = 1.54;  % Units [µm, microns]
lambda_max = 1.56;
lambda_step = 0.01e-3; % wavelength step [microns]
                       % Typical minimum step for a tunable laser is 1-10 pm.
lambda=lambda_min:lambda_step:lambda_max;


% Define the waveguide effective index compact model:
% - as a Taylor expansion around the central wavelength, lambda0
%  use a Matlab anonymous function 
% lambda0 = 1.55; a1=2.4489; a2=-1.1337; a3=-0.0451;  % these are constants from the waveguide model.
lambda0 = 1.55; a1=2.4; a2=0; a3=-0.0;  % these are constants from the waveguide model.

neff = @(lambda) ...
		(a1 + a2.*(lambda-lambda0) + a3.*(lambda-lambda0).^2); 
        
% plot the effective index version wavelength, and check if this is as expected:
figure;
plot(lambda, neff(lambda),'LineWidth',3);

% Complex propagation constant for the waveguide
% alpha = 1e-3;  % propagation loss [micron^-1]; constant
alpha = 0e-3;  % propagation loss [micron^-1]; constant
beta = @(lambda) ...
		(2*pi*neff(lambda)./lambda - 1i*alpha/2*ones(1,length(lambda)) );


% Define the MZI transfer function
T_MZI = @(L1, L2, lambda) ...
        ( 0.25* abs(exp(-1i*beta(lambda)*L1)+exp(-1i*beta(lambda)*L2)).^2);

% Define the two waveguide lengths in the MZI
L1=100;  % Waveguide 1 length, Units [µm, microns]
L2=100.0953;  % Waveguide 2 length, Units [µm, microns]

% plot the MZI transfer function, and check if this is as expected:
% plot in linear scale:
figure;
plot(lambda, T_MZI(L1, L2, lambda),'LineWidth',3);
xlabel ('Wavelength [\mum]');
ylabel ('Transmission');
axis tight
title ('MZI transfer function');

% plot in log (dB) scale:
figure;
T_MZI_dB = 10*log10(T_MZI(L1, L2, lambda));
plot(lambda, T_MZI_dB,'LineWidth',3);
xlabel ('Wavelength [\mum]');
ylabel ('Transmission [dB]');
axis tight
title ('MZI transfer function');
      