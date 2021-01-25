% Phot1x_MZI_matlab.m
% by Lukas Chrostowski, 2015

% the wavelength range of interest.
lambda_min = 1.5;  % Units [µm, microns]
lambda_max = 1.6;
lambda_step = 0.01e-3; % wavelength step [microns]
                       % Typical minimum step for a tunable laser is 1-10 pm.
lambda=lambda_min:lambda_step:lambda_max;

% Define the MZI transfer function
%  use Matlab anonymous functions

% Effective index:
% - as a Taylor expansion around the central wavelength, lambda0
lambda0 = 1.55; n1=2.4489; n2=--1.1337; n3=-0.0451;  % these are constants from the waveguide model.
neff = @(lambda) ...
		(n1 + n2.*(lambda-lambda0) + n3.*(lambda-lambda0).^2); 
% plot, and check if this is as expected:
figure;
plot(lambda, neff(lambda),'LineWidth',3);

% Complex propagation constant
alpha_dBpercm = 5;
alpha = alpha_dBpercm / 4.34 * 1e-4;  % propagation loss [micron^-1]; constant
beta = @(lambda) ...
		(2*pi*neff(lambda)./lambda - 1i*alpha/2*ones(1,length(lambda)) );


% MZI transfer function
T_MZI = @(L1, L2, lambda) ...
        ( 0.25* abs(exp(-1i*beta(lambda)*L1)+exp(-1i*beta(lambda)*L2)).^2);

% plot, and check if this is as expected:
L1=100;
L2=200;  % Units [µm, microns], variable
figure;
plot(lambda, T_MZI(L1, L2, lambda),'LineWidth',3);
xlabel ('Wavelength [\mum]');
ylabel ('Transmission');
axis tight
title ('MZI transfer function');

figure;
T_MZI_dB = 10*log10(T_MZI(L1, L2, lambda));
plot(lambda, T_MZI_dB,'LineWidth',3);
xlabel ('Wavelength [\mum]');
ylabel ('Transmission [dB]');
axis tight
title ('MZI transfer function');
