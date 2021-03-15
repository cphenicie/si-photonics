simulated_n = [2.4489, -1.1337, -0.0451];

corner_file_dir = 'C:\Users\chris\Documents\GitHub\si-photonics\phot1x\05_fabrication\corner_analysis';
corner_file_names = {%'wg_2D_sweep_wl_width=470nm,thick=215.3nm.mat',
                     'wg_2D_sweep_wl_width=470nm,thick=223.1nm.mat',
                     'wg_2D_sweep_wl_width=510nm,thick=215.3nm.mat'
                     %'wg_2D_sweep_wl_width=510nm,thick=223.1nm.mat'
                     };
                             
corner_legend_names = {%'(470nm x 215.3nm)',
                       '(470nm x 223.1nm)',
                       '(510nm x 215.3nm)'
                       %'(510nm x 223.1nm)'
                       };
                   
fit_file_dir = 'C:\Users\chris\Documents\GitHub\si-photonics\phot1x\07_final_analysis-exam';
fit_file_names = {'MZI30um_fit.mat', 
                  'MZI50um_fit.mat', 
                  'MZI100um_fit.mat',
                  'MZI300um_fit.mat',
                  'MZI500um_fit.mat',
                  'MIA_fit.mat',
                  'MIB_fit.mat'
                  };

fit_legend_names = {'MZI 30um', 
                    'MZI 50um', 
                    'MZI 100um', 
                    'MZI 300um', 
                    'MZI 500um',
                    'MIA',
                    'MIB'
                    };
                
mi_file_names = {'MIA_fit.mat', 
                  'MIB_fit.mat'
                  };

mi_legend_names = {'MI A', 
                    'MI B' 
                    };


n1_vals = zeros(1,length(fit_file_names));
n2_vals = zeros(1,length(fit_file_names));
n3_vals = zeros(1,length(fit_file_names));


n1_errs = zeros(1,length(fit_file_names));
n2_errs = zeros(1,length(fit_file_names));
n3_errs = zeros(1,length(fit_file_names));

ng_vals = zeros(1,length(fit_file_names));
ng_errs = zeros(1,length(fit_file_names));

fsr_vals = zeros(1,length(fit_file_names));
fsr_errs = zeros(1,length(fit_file_names));

%% Plot corner analysis and data for MZIs

figure();
hold on

lambda0 = 1550e-9;
lambda_fit = 1530e-9:1e-10:1570e-9;

neff_sim = neff(lambda0*1e6, simulated_n, lambda_fit*1e6);
ng_sim = ng_func(neff_sim, lambda_fit*1e6);
plot(lambda_fit*1e6, ng_sim, 'k--', 'DisplayName', 'Simulated waveguide', 'Linewidth', 2)


for i=1:length(corner_file_names)

   load(fullfile(corner_file_dir, corner_file_names{i}))
   n1 = real(a_fit(1));
   n2 = real(a_fit(2));
   n3 = real(a_fit(3));
 
   neff_corner = neff(lambda0*1e6, [n1, n2, n3], lambda_fit*1e6);
   ng_corner = ng_func(neff_corner, lambda_fit*1e6);
   label = join(['Corner ', corner_legend_names{i}]);
   plot(lambda_fit*1e6, ng_corner, '--', 'DisplayName', label, 'Linewidth', 2)
end

for i=1:length(fit_file_names)

   load(fullfile(fit_file_dir, fit_file_names{i}))
   n1 = xfit(1);
   n2 = xfit(2);
   n3 = xfit(3);
 
   n1_vals(i) = n1;
   n2_vals(i) = n2;
   n3_vals(i) = n3;
   
   n1_errs(i) = xerr(1);
   n2_errs(i) = xerr(2);
   n3_errs(i) = xerr(3);
   
   neff_fit = neff(lambda0*1e6, [n1, n2, n3], lambda_fit*1e6);
   ng_fit = ng_func(neff_fit, lambda_fit*1e6);
   label = join(['Fit ', fit_legend_names{i}]);
   plot(lambda_fit*1e6, ng_fit, 'DisplayName', label, 'Linewidth', 1.5)
   
   ng_vals(i) = ng_fit(lambda_fit*1e6 == 1.55);
   ng_errs(i) = sqrt(xerr(1)^2 + xerr(2)^2*lambda0*1e6);
   
   fsr_vals(i) = fsr;
   fsr_errs(i) = fsr_err;
   
   
end

legend();
xlabel('Wavelength [\mum]')
ylabel('Group index')

xlim([min(lambda_fit)*1e6, max(lambda_fit)*1e6])

grid on
box on

format long

[n1_vals', n2_vals', n3_vals', ng_vals', 1000*fsr_vals']

[n1_errs', n2_errs', n3_errs', ng_errs', 1000*fsr_errs']

n1_weights = 1./n1_errs.^2;
n2_weights = 1./n2_errs.^2;
n3_weights = 1./n3_errs.^2;

% n1_avg = sum(n1_vals .* n1_weights) / sum(n1_weights);
% n2_avg = sum(n2_vals .* n2_weights) / sum(n2_weights);
% n3_avg = sum(n3_vals .* n3_weights) / sum(n3_weights);
n1_avg = mean(n1_vals);
n2_avg = mean(n2_vals);
n3_avg = mean(n3_vals);

% n1_std = 1/sqrt(sum(n1_weights));
% n2_std = 1/sqrt(sum(n2_weights));
% n3_std = 1/sqrt(sum(n3_weights));
n1_std = std(n1_vals);
n2_std = std(n2_vals);
n3_std = std(n3_vals);

sprintf('n1_vals = %.3f +/- %.3f', n1_avg, n1_std)
sprintf('n2_vals = %.3f +/- %.3f', n2_avg, n2_std)
sprintf('n3_vals = %.3f +/- %.3f', n3_avg, n3_std)
sprintf('ng_vals = %.3f +/- %.3f', mean(ng_vals), std(ng_vals))

%% Plot corner analysis and data for MIs

figure();
hold on

lambda0 = 1550e-9;
lambda_fit = 1530e-9:1e-10:1570e-9;

neff_sim = neff(lambda0*1e6, simulated_n, lambda_fit*1e6);
ng_sim = ng_func(neff_sim, lambda_fit*1e6);
plot(lambda_fit*1e6, ng_sim, 'k--', 'DisplayName', 'Simulated waveguide', 'Linewidth', 2)


for i=1:length(corner_file_names)

   load(fullfile(corner_file_dir, corner_file_names{i}))
   n1 = real(a_fit(1));
   n2 = real(a_fit(2));
   n3 = real(a_fit(3));
 
   neff_corner = neff(lambda0*1e6, [n1, n2, n3], lambda_fit*1e6);
   ng_corner = ng_func(neff_corner, lambda_fit*1e6);
   label = join(['Corner ', corner_legend_names{i}]);
   plot(lambda_fit*1e6, ng_corner, '--', 'DisplayName', label, 'Linewidth', 2)
end

for i=1:length(mi_file_names)

   load(fullfile(fit_file_dir, mi_file_names{i}))
   n1 = xfit(1);
   n2 = xfit(2);
   n3 = xfit(3);
 
   n1_vals(i) = n1;
   n2_vals(i) = n2;
   n3_vals(i) = n3;
   
   n1_errs(i) = xerr(1);
   n2_errs(i) = xerr(2);
   n3_errs(i) = xerr(3);
   
   neff_fit = neff(lambda0*1e6, [n1, n2, n3], lambda_fit*1e6);
   ng_fit = ng_func(neff_fit, lambda_fit*1e6);
   label = join(['Fit ', mi_legend_names{i}]);
   plot(lambda_fit*1e6, ng_fit, 'DisplayName', label, 'Linewidth', 1.5)
   
   ng_vals(i) = ng_fit(lambda_fit*1e6 == 1.55);
   ng_errs(i) = sqrt(xerr(1)^2 + xerr(2)^2*lambda0*1e6);
   
   join([mi_legend_names{i}, ' fsr=', num2str(fsr*1e6), ' pm, dL = ',...
       num2str(dL), ' um, expected ng from FSR alone: ',... 
       num2str(1.55^2 / (fsr * dL))])
   
   
end

legend()
xlabel('Wavelength [\mum]')
ylabel('Group index')

xlim([min(lambda_fit)*1e6, max(lambda_fit)*1e6])

grid on
box on

%% Some useful functions

% effective index:
function n = neff(lambda0, nx, lambda)
    n = (nx(1) + nx(2).*(lambda-lambda0) + nx(3).*(lambda-lambda0).^2);
end

% Group index, from finite differences of the effective index
function ng_arr = ng_func(neff, lambda)
    % lambda: array of wavelengths
    % neff: array of effective indices evaluated at each wavelength
    dndlambda = diff(neff) ./ diff(lambda); 
    dndlambda = [dndlambda, dndlambda(end)]; % Make derivative same size
    ng_arr = (neff - lambda .* dndlambda);
end