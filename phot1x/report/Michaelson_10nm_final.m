load C:/Users/chris/Documents/GitHub/si-photonics/phot1x/report/Michaelson_10nm_final.mat;
figure(1);
ax1 = axes; 
plot(lum.x0, lum.y0,lum.x1, lum.y1)
set(ax1, 'XLim', [1548.5 1548.8])
set(ax1, 'YLim', [0 0.2])
set(ax1,'XGrid', 'on')
set(ax1,'YGrid', 'on')
legend('width = 500 nm', 'width = 510 nm')
xlabel('Wavelength [nm]')
ylabel('Transmission [I_o / I_i]')

figure(2);
ax1 = axes; 
plot(lum.x0, lum.y0,lum.x1, lum.y1)
set(ax1, 'XLim', [1552.45 1552.75])
set(ax1, 'YLim', [0 0.2])
set(ax1,'XGrid', 'on')
set(ax1,'YGrid', 'on')
legend('width = 500 nm', 'width = 510 nm')
xlabel('Wavelength [nm]')
ylabel('Transmission [I_o / I_i]')

% User inputs
delta_L = 2 * 10e6; % nm, factor of 2 because Michaelson interferometer
first_wl = 1540; % nm
last_wl = 1560; % nm
data_x_all_1 = lum.x0;
data_y_all_1 = lum.y0;
data_x_all_2 = lum.x1;
data_y_all_2 = lum.y1;



good_idx = find(data_x_all_1 > first_wl & data_x_all_1 < last_wl);
data_x_1 = data_x_all_1(good_idx);
data_y_1 = data_y_all_1(good_idx);
[pks, locs] = findpeaks(data_y_1);
max_wl = data_x_1(locs);
lambda = mean(data_x_1);
fsr_1 = abs(mean(diff(max_wl)));
ng_1 = lambda^2 / (fsr_1 * delta_L);

good_idx = find(data_x_all_2 > first_wl & data_x_all_2 < last_wl);
data_x_2 = data_x_all_2(good_idx);
data_y_2 = data_y_all_2(good_idx);
[pks, locs] = findpeaks(data_y_2);
max_wl = data_x_2(locs);
lambda = mean(data_x_2);
fsr_2 = abs(mean(diff(max_wl)));
ng_2 = lambda^2 / (fsr_2 * delta_L);

delta_fsr = fsr_1 - fsr_2;
dngdw = -1.69e-3; % Computed separately
delta_w = delta_fsr * delta_L * ng_1^2 / (lambda^2 * dngdw)