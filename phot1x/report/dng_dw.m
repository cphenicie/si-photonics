c = 299792458.0; % Speed of light in m/s
load C:/Users/chris/Documents/GitHub/si-photonics/phot1x/report/490um_waveguide.mat;
lambda = c ./ f * 1e9; 
ng_490 = c ./ vg;

load C:/Users/chris/Documents/GitHub/si-photonics/phot1x/report/495um_waveguide.mat;
ng_495 = c ./ vg;

load C:/Users/chris/Documents/GitHub/si-photonics/phot1x/report/500um_waveguide.mat;
ng_500 = c ./ vg;

load C:/Users/chris/Documents/GitHub/si-photonics/phot1x/report/505um_waveguide.mat;
ng_505 = c ./ vg;

load C:/Users/chris/Documents/GitHub/si-photonics/phot1x/report/510um_waveguide.mat;
ng_510 = c ./ vg;

figure(1)
plot(lambda, ng_490, lambda, ng_495, lambda, ng_500, lambda, ng_505, lambda, ng_510)
xlabel('Wavelength [nm]')
ylabel('Group index')
legend('490nm', '495nm', '500nm', '505nm', '510nm')

width = [490, 495, 500, 505, 510];
idx = 6;
ng = [ng_490(idx), ng_495(idx), ng_500(idx), ng_505(idx),ng_510(idx)];

figure(2)
plot(width, ng, 'o')

width_0 = 500;
ng_eq = @(nx, width) ...
		(nx(1) + nx(2).*(width-width_0)); 
X=[4.2, -0.01]; 

% curve fit to find expression for ng.
format long
X = lsqcurvefit (ng_eq, X, width, ng)

hold on 
plot(width, ng_eq(X, width))

legend('data', 'linear fit')
xlabel('Waveguide width [nm]')
ylabel('Group index at 1554nm')
title(sprintf('dng/dw = %.3e', X(2)))
