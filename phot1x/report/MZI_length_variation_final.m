load C:/Users/chris/Documents/GitHub/si-photonics/phot1x/report/MZI_length_variation_final.mat;
ax1 = axes; 
% plot(lum.x0, lum.y0,lum.x1, lum.y1,lum.x2, lum.y2,lum.x3, lum.y3,...
%    lum.x4,lum.y4, 'Linewidth', 1)
plot(lum.x0, lum.y0, lum.x4,lum.y4, 'Linewidth', 1)
set(ax1, 'XLim', [1524 1576])
set(ax1, 'YLim', [0 0.35])
set(ax1,'XGrid', 'on')
set(ax1,'YGrid', 'on')
%legend('\DeltaL = 30 \mum','\DeltaL = 50 \mum','\DeltaL = 100 \mum',...
%'\DeltaL = 300 \mum','\DeltaL = 500 \mum')
legend('\DeltaL = 30 \mum', '\DeltaL = 500 \mum')
xlabel('Wavelength [nm]')
ylabel('Transmission [I_o / I_i]')

lambda = 1550;
n_g = 4.206;
delta_L_list = [30e3, 50e3, 100e3, 300e3, 500e3];
fsr_list = 1550^2 ./ (delta_L_list * n_g)


% User inputs
delta_L = 500e3;
first_wl = 1500;
last_wl = 1600;

x_list = {lum.x0, lum.x1, lum.x2, lum.x3, lum.x4};
y_list = {lum.y0, lum.y1, lum.y2, lum.y3, lum.y4};

for i = 1:5
    data_x_all = x_list{i};
    data_y_all = y_list{i};
    good_idx = find(data_x_all > first_wl & data_x_all < last_wl);
    data_x = data_x_all(good_idx);
    data_y = data_y_all(good_idx);
    % plot(data_x, data_y)
    % hold on

    [pks, locs] = findpeaks(data_y);
    max_wl = data_x(locs);
    % plot(max_wl, pks, 'o')

    lambda = mean(data_x);
    fsr = abs(mean(diff(max_wl)))
%     ng = lambda^2 / (fsr * delta_L)
    ng = lambda^2 / (fsr * delta_L_list(i))

%     xlabel('Wavelength [nm]')
%     ylabel('Transmission [I_o / I_i]')
%     title(sprintf('Given dL = %.1f um, Measured FSR = %.2f nm, Computed ng = %.3f',...
%         delta_L/1000, fsr, ng))
% 
%     legend('Sim. Data', 'Extracted peaks')
end