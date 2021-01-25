load C:/Users/chris/Documents/GitHub/si-photonics/phot1x/report/MZI_length_variation.mat;
ax1 = axes; 
plot(lum.x0, lum.y0, 'k', lum.x1, lum.y1, 'r', lum.x2, lum.y2, 'b', 'Linewidth', 1)
set(ax1, 'XLim', [1500 1600])
set(ax1, 'YLim', [0 1])
set(ax1,'XGrid', 'on')
set(ax1,'YGrid', 'on')
legend('\DeltaL = 5 \mum (1.1x)','\DeltaL = 50 \mum (2x)','\DeltaL = 200 \mum (5x)')
xlabel('Wavelength [nm]')
ylabel('Transmission [I_o / I_i]')
