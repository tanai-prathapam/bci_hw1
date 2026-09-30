[eeg,header] = sload('Subject_006_TESS_Online__feedback__s006_r001_2021_08_30_163523.gdf');

%% temporal filtering
[b,a] = butter(order,norm_cutoff_freq,type)

norm_cutoff_freq = [low high]/Nyquist 

filt(b,a,data) or 
filtfilt(b,a,data)

%% spatial filtering

CAR or Laplacian

%% trial extraction

data [time samples x channels x trials]

%% mu power 

instantenous power = time samples ^ 2

for each trial
    average power = sum of instanteneous power/n of samples

%% topoplot
topoplot(data vector,selectedChannels) % data vector [n of channels x 1]

%% grand average plot

%% CCA
% maximize correlation between the single-trial data and the class grand average template

[spatialFilter,~] =  canoncorr(concat_data', concat_ga'); % spatial filter [n channels x n components]