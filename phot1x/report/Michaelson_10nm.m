load C:/Users/chris/Documents/GitHub/si-photonics/phot1x/report/Michaelson_10nm.mat;
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
set(ax1, 'XLim', [1552.3 1552.6])
set(ax1, 'YLim', [0 0.2])
set(ax1,'XGrid', 'on')
set(ax1,'YGrid', 'on')
legend('width = 500 nm', 'width = 510 nm')
xlabel('Wavelength [nm]')
ylabel('Transmission [I_o / I_i]')