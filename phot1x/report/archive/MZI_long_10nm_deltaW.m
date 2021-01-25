load C:/Users/chris/Documents/GitHub/si-photonics/phot1x/report/MZI_long_10nm_deltaW.mat;
ax1 = axes; 
plot(lum.x0, lum.y0,lum.x1, lum.y1)
set(ax1, 'XLim', [1524.53 1576.78])
set(ax1, 'YLim', [-0.0690662 0.27852])
set(ax1,'XGrid', 'on')
set(ax1,'YGrid', 'on')
xlabel('Wavelength [nm]')
ylabel('Transmission [I_o / I_i]')
legend('w = 500nm', 'w = 510 nm')
