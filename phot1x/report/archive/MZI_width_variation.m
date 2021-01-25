load C:/Users/chris/Documents/GitHub/si-photonics/phot1x/report/MZI_width_variation.mat;
ax1 = axes; 
plot(lum.x1, lum.y1, 'k', lum.x2, lum.y2, 'r', lum.x0, lum.y0, 'b', 'Linewidth', 1)
set(ax1, 'XLim', [1500 1600])
set(ax1, 'YLim', [0 1])
set(ax1,'XGrid', 'on')
set(ax1,'YGrid', 'on')
legend('\Deltaw = 0.05 \mum (1.1x)', '\Deltaw = 0.5 \mum (2x)', '\Deltaw = 2 \mum (5x)')
xlabel('Wavelength [nm]')
ylabel('Transmission [I_o / I_i]')
