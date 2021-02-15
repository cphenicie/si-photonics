file_dir = 'C:\Users\chris\Documents\GitHub\si-photonics\phot1x\05_fabrication\corner_analysis';
file_names = {'wg_2D_sweep_wl_width=470nm,thick=215.3nm.mat',
              'wg_2D_sweep_wl_width=470nm,thick=223.1nm.mat',
              'wg_2D_sweep_wl_width=510nm,thick=215.3nm.mat',
              'wg_2D_sweep_wl_width=510nm,thick=223.1nm.mat',
              'wg_2D_sweep_wl_width=500nm,thick=220nm.mat'};
          
width = [470, 470, 510, 510, 500];
thick = [215.3, 223.1, 215.3, 223.1, 220];

legend_names = {'width=470nm,thick=215.3nm',
                'width=470nm,thick=223.1nm',
                'width=510nm,thick=215.3nm',
                'width=510nm,thick=223.1nm',
                'width=510nm,thick=215.3nm'};

L1=100;
L2=200;  % Units [µm], variable

lambda0 = 1.55;      
lambda_min = 1.54;  % Units [Âµm, microns]
lambda_max = 1.56;
lambda_step = 0.01e-3; % wavelength step [microns]
                       % Typical minimum step for a tunable laser is 1-10 pm.
lambda=lambda_min:lambda_step:lambda_max;
           

                
figure;
hold on
n_files = length(file_names);
for i=1:n_files

   load(fullfile(file_dir, file_names{i}))
   n1 = real(a_fit(1));
   n2 = real(a_fit(2));
   n3 = real(a_fit(3));
 
   neff_cm = @(lambda) ...
		(n1 + n2.*(lambda-lambda0) + n3.*(lambda-lambda0).^2); 
    
   alpha_dBpercm = 5;
   alpha = alpha_dBpercm / 4.34 * 1e-4;  % propagation loss [micron^-1]; constant
   beta = @(lambda) ...
           (2*pi*neff_cm(lambda)./lambda - 1i*alpha/2*ones(1,length(lambda)));
   T_MZI = @(L1, L2, lambda) ...
           (0.25*abs(exp(-1i*beta(lambda)*L1) + exp(-1i*beta(lambda)*L2)).^2);

   data_y = T_MZI(L1, L2, lambda);
   plot(lambda, 10*log10(data_y),'LineWidth',0.1);
   
   [pks, locs] = findpeaks(data_y);
   max_wl = lambda(locs);
   % plot(max_wl, pks, 'o')

   fsr = abs(mean(diff(max_wl))) * 1e3;
   sprintf(join(['width=', num2str(width(i)), 'nm',  ...
        '; thickness=', num2str(thick(i)), 'nm', ...
        '; fsr = ', num2str(fsr), 'nm', ...
        '; neff = ', num2str(neff_cm(lambda0)),...
        '; ng = ', num2str(ng_mode)]))

end

xlabel ('Wavelength [\mum]');
ylabel ('Transmission [dB]');
axis tight
title ('MZI transfer function');
legend(legend_names)