
% Read data files from experiments
% Enter the Dropbox URLs here.  Make sure the URL has a =1 at the end:
%  Loopback structure:
    loopback_name = 'Phot1x-DataSets/DataSet1/ZiheGao_MZI1_272_Scan1.mat';
% 	url_loopback = 'https://www.dropbox.com/s/w915qfix9kwlwv7/ZiheGao_MZI1_272_Scan1.mat?dl=1';
%  MZI:
    mzi_name = 'Phot1x-DataSets/DataSet1/ZiheGao_MZI2_271_Scan1.mat';
% 	url_mzi = 'https://www.dropbox.com/s/1rvjfef4jqybc12/ZiheGao_MZI2_271_Scan1.mat?dl=1';
% Calibrate the MZI data using the loopback structure
% Plot


PORT=1; % Which Fibre array port is the output connected to?
FONTSIZE=20;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loopback data:
% a=websave('loopback.mat',url_loopback); % get data from Dropbox
% load('loopback.mat');
load(loopback_name);
% Data is stored in variable "scanResults".
% There are two columns - wavelength (1), and amplitude (2)
lambda=scanResults(1,PORT).Data(:,1)/1e9;
amplitude=scanResults(1,PORT).Data(:,2);
figure;
plot (lambda*1e6, amplitude);
title ('Calibration loopback'); 
xlabel ('Wavelength [\mum]','FontSize',FONTSIZE)
ylabel ('Insertion Loss [dB]','FontSize',FONTSIZE)
hold all;

% Fit the data with a polynomial
p=polyfit((lambda-mean(lambda))*1e6, amplitude, 5);
amplitude_LOOPBACK=polyval(p,(lambda-mean(lambda))*1e6);
plot (lambda*1e6, amplitude_LOOPBACK);
% find wavelength range with usable data, in the loopback
loopback_IL = max(amplitude);
new_lambda_i=find(amplitude>loopback_IL-10);
lambda=lambda(new_lambda_i);
lambda_min = min(lambda);
lambda_max = max(lambda);
amplitude=amplitude(new_lambda_i);
% refit the loopback
LOOPBACK=polyfit((lambda-mean(lambda))*1e6, amplitude, 4);
amplitude_LOOPBACK=polyval(LOOPBACK,(lambda-mean(lambda))*1e6);
plot (lambda*1e6, [amplitude_LOOPBACK],'r-','Linewidth',1);
axis tight;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MZI data:
% a=websave('mzi.mat',url_mzi); % get data from Dropbox
% load('mzi.mat');
load(mzi_name);
lambda1=scanResults(1,PORT).Data(:,1)/1e9;
amplitude=scanResults(1,PORT).Data(:,2);
figure;
plot (lambda1*1e6, amplitude);
% title ('MZI (raw data)'); 
% xlabel ('Wavelength [\mum]','FontSize',FONTSIZE)
% ylabel ('Insertion Loss [dB]','FontSize',FONTSIZE)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MZI data - calibrated
%
% data only within the bandwidth of interest.
lambda=lambda_min:min(diff(lambda1)):lambda_max;
amplitude=interp1(lambda1, amplitude, lambda,'linear');
amplitude(find(amplitude==-inf))=-50;
% calibrate data
amplitude_cal=amplitude-polyval(LOOPBACK,(lambda-mean(lambda))*1e6);
hold on;
% figure;
plot (lambda*1e6, amplitude_cal);
title ('MZI (calibrated with loopback)'); 
xlabel ('Wavelength [\mum]','FontSize',FONTSIZE)
ylabel ('Insertion Loss [dB]','FontSize',FONTSIZE)
legend('Raw', 'Calibrated')

    
      