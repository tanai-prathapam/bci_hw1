%% 0. Setup
clear; clc; close all;
biosigDir = "C:\Users\bravo\OneDrive - The University of Texas at Austin\Desktop\Tanai's Data\UT Undergrad ALL\UT Classes ALL\BCI\HW_1_export\toolboxes\toolboxes\biosig";
assert(isfolder(biosigDir), 'BioSig folder not found');
addpath(genpath(biosigDir));

% EDIT if your data folder is somewhere else
dataDir = fullfile(pwd, 'Subject_006_Session_006_TESS_Online_Visual');
assert(isfolder(dataDir), 'Data folder not found: %s', dataDir);

files = dir(fullfile(dataDir, '*.gdf'));
[~, ord] = sort({files.name});  files = files(ord);     % r001..r004 in order
nRuns = numel(files);
fprintf('\n=== 0. Setup ===\n');
fprintf('Found %d GDF files in %s\n', nRuns, dataDir);
for r = 1:nRuns, fprintf('  run %d: %s\n', r, files(r).name); end
assert(nRuns == 4, 'Expected 4 runs, found %d', nRuns);

%% 1. Parameters (justify the choices in your report)
order = 4;  band = [8 12];                 % mu band, Hz
keepCh = 1:32;                             % last 4 columns ignored
% Event codes: VERIFY against biosig\eventcodes.txt before trusting
codeLH = hex2dec('0301');                  % assumed left-hand cue
codeRH = hex2dec('0302');                  % assumed right-hand cue
taskDur_s = NaN;                           % task duration after cue, s (set after inspecting events)
lastWin_s = 0.5;                           % last 0.5 s of task for topoplots

%% 2. Channel labels, Laplacian neighbors, spatial filter matrices
chanLabels = {'FP1','FPZ','FP2','F7','F3','FZ','F4','F8','FC5','FC1','FC2','FC6', ...
              'M1','T7','C3','CZ','C4','T8','M2','CP5','CP1','CP2','CP6','P7','P3', ...
              'PZ','P4','P8','POZ','O1','OZ','O2'};
assert(numel(chanLabels) == 32);

nbrNames = struct( ...
 'FP1', {{'FPZ','F3','F7'}}, ...
 'FPZ', {{'FP1','FP2','FZ'}}, ...
 'FP2', {{'FPZ','F4','F8'}}, ...
 'F7',  {{'FP1','F3','FC5'}}, ...
 'F3',  {{'FP1','F7','FZ','FC5','FC1'}}, ...
 'FZ',  {{'FPZ','F3','F4','FC1','FC2'}}, ...
 'F4',  {{'FP2','FZ','F8','FC2','FC6'}}, ...
 'F8',  {{'FP2','F4','FC6'}}, ...
 'FC5', {{'F7','F3','FC1','T7','C3'}}, ...
 'FC1', {{'F3','FZ','FC5','FC2','C3','CZ'}}, ...
 'FC2', {{'FZ','F4','FC1','FC6','CZ','C4'}}, ...
 'FC6', {{'F4','F8','FC2','C4','T8'}}, ...
 'M1',  {{'T7','CP5'}}, ...
 'T7',  {{'FC5','C3','CP5'}}, ...
 'C3',  {{'FC5','FC1','T7','CZ','CP5','CP1'}}, ...
 'CZ',  {{'FC1','FC2','C3','C4','CP1','CP2'}}, ...
 'C4',  {{'FC2','FC6','CZ','T8','CP2','CP6'}}, ...
 'T8',  {{'FC6','C4','CP6'}}, ...
 'M2',  {{'T8','CP6'}}, ...
 'CP5', {{'T7','C3','CP1','P7','P3'}}, ...
 'CP1', {{'C3','CZ','CP5','CP2','P3','PZ'}}, ...
 'CP2', {{'CZ','C4','CP1','CP6','PZ','P4'}}, ...
 'CP6', {{'C4','T8','CP2','P4','P8'}}, ...
 'P7',  {{'CP5','P3','O1'}}, ...
 'P3',  {{'CP5','CP1','P7','PZ','O1','POZ'}}, ...
 'PZ',  {{'CP1','CP2','P3','P4','POZ'}}, ...
 'P4',  {{'CP2','CP6','PZ','P8','POZ','O2'}}, ...
 'P8',  {{'CP6','P4','O2'}}, ...
 'POZ', {{'PZ','P3','P4','OZ'}}, ...
 'O1',  {{'P7','P3','POZ','OZ'}}, ...
 'OZ',  {{'POZ','O1','O2'}}, ...
 'O2',  {{'P8','P4','POZ','OZ'}});

nbr = cell(32,1);
for i = 1:32
    names = nbrNames.(upper(chanLabels{i}));
    nbr{i} = cellfun(@(n) find(strcmpi(chanLabels, n)), names);
end
assert(all(cellfun(@numel, nbr) >= 1), 'A channel has no neighbors');
assert(~any(arrayfun(@(i) ismember(i, nbr{i}), 1:32)), 'A channel lists itself');

W_car = eye(32) - ones(32)/32;
W_lap = eye(32);
for i = 1:32
    W_lap(i, nbr{i}) = -1/numel(nbr{i});
end
fprintf('\n=== 2. Spatial filters defined ===\n');
fprintf('Neighbors per channel:\n');
for i = 1:32
    fprintf('  %-4s <- %s\n', chanLabels{i}, strjoin(chanLabels(nbr{i}), ', '));
end
%% 3. Load and filter every run
EEG = cell(nRuns,1);  HDR = cell(nRuns,1);
MU  = cell(nRuns,1);  CAR = cell(nRuns,1);  LAP = cell(nRuns,1);
for r = 1:nRuns
    fn = fullfile(dataDir, files(r).name);
    [s, h] = sload(fn);
    fs_r = h.SampleRate;
    if r == 1, fs = fs_r; chanLabels = h.Label(keepCh); end

    assert(fs_r == fs, 'Sampling rate differs in run %d', r);
    assert(isequal(h.Label(keepCh), chanLabels), 'Channel labels differ in run %d', r);
    assert(size(s,2) >= 32, 'Run %d has fewer than 32 columns', r);

    eeg = s(:, keepCh);
    assert(~any(isnan(eeg(:))), 'NaN in raw data, run %d', r);

    [b, a] = butter(order, band/(fs/2), 'bandpass');
    mu = filtfilt(b, a, eeg);

    EEG{r} = eeg;  HDR{r} = h;
    MU{r}  = mu;
    CAR{r} = mu * W_car';
    LAP{r} = mu * W_lap';

    fprintf('\n--- Run %d: %s ---\n', r, files(r).name);
    fprintf('Raw size %d x %d -> kept %d x %d | fs %g Hz | length %.1f s\n', ...
        size(s,1), size(s,2), size(eeg,1), size(eeg,2), fs_r, size(eeg,1)/fs_r);
    fprintf('Dropped columns: %s\n', strjoin(h.Label(33:end), ', '));
    fprintf('Filter stable: %d | NaN/Inf after filter: %d\n', ...
        all(abs(roots(a)) < 1), any(~isfinite(mu(:))));
    fprintf('Ch1 (%s) mu-band power fraction: before %.3f, after %.3f\n', chanLabels{1}, ...
        bandpower(eeg(:,1), fs, band)/bandpower(eeg(:,1), fs, [0.5 fs/2-1]), ...
        bandpower(mu(:,1),  fs, band)/bandpower(mu(:,1),  fs, [0.5 fs/2-1]));
    fprintf('CAR max |mean across ch|: %.2e | Laplacian NaN: %d\n', ...
        max(abs(mean(CAR{r},2))), any(isnan(LAP{r}(:))));
end

% Check that the matrix form of the Laplacian equals the explicit loop (run 1)
lap_loop = MU{1};
for i = 1:32, lap_loop(:,i) = MU{1}(:,i) - mean(MU{1}(:, nbr{i}), 2); end
fprintf('\nMatrix vs loop Laplacian max diff (run 1): %.2e\n', max(abs(lap_loop(:) - LAP{1}(:))));

%% 4. Inspect events in every run
fprintf('\n=== 4. Events per run ===\n');
for r = 1:nRuns
    typ = HDR{r}.EVENT.TYP;  pos = HDR{r}.EVENT.POS;
    [u,~,k] = unique(typ);  cnt = accumarray(k,1);
    fprintf('Run %d: %d events\n', r, numel(typ));
    for i = 1:numel(u)
        fprintf('   code 0x%04X (%d): %d events\n', u(i), u(i), cnt(i));
    end
    n = min(12, numel(typ));
    fprintf('   First %d events (code, time s):\n', n);
    disp([typ(1:n), pos(1:n)/fs]);
    if isfield(HDR{r}.EVENT,'DUR') && ~isempty(HDR{r}.EVENT.DUR)
        fprintf('   First durations (samples): %s\n', mat2str(HDR{r}.EVENT.DUR(1:n)'));
    end
end

%% 5. Count trials per class and run
fprintf('\n=== 5. Trial counts ===\n');
nLH = zeros(nRuns,1);  nRH = zeros(nRuns,1);
for r = 1:nRuns
    typ = HDR{r}.EVENT.TYP;
    nLH(r) = sum(typ == codeLH);  nRH(r) = sum(typ == codeRH);
    fprintf('Run %d: LH %d, RH %d\n', r, nLH(r), nRH(r));
end
fprintf('Total: LH %d, RH %d, all %d\n', sum(nLH), sum(nRH), sum(nLH)+sum(nRH));
if sum(nLH)+sum(nRH) == 0
    warning('No trials found with the assumed codes. Fix codeLH/codeRH from eventcodes.txt.');
end

%% 6. Extract trials and compute mu power
if isnan(taskDur_s) || sum(nLH)+sum(nRH) == 0
    fprintf('\nSet taskDur_s and valid class codes above, then rerun sections 6 onward.\n');
else
    nS = round(taskDur_s*fs);
    cls = {'LH','RH'};  codes = [codeLH codeRH];
    T = struct();
    for c = 1:2
        T.(cls{c}).mu = [];  T.(cls{c}).car = [];  T.(cls{c}).lap = [];  T.(cls{c}).run = [];
        for r = 1:nRuns
            typ = HDR{r}.EVENT.TYP;  pos = HDR{r}.EVENT.POS;
            st = pos(typ == codes(c));
            for t = 1:numel(st)
                idx = st(t) : st(t)+nS-1;
                if idx(end) > size(MU{r},1), continue; end    % skip truncated trial
                T.(cls{c}).mu  = cat(3, T.(cls{c}).mu,  MU{r}(idx,:));
                T.(cls{c}).car = cat(3, T.(cls{c}).car, CAR{r}(idx,:));
                T.(cls{c}).lap = cat(3, T.(cls{c}).lap, LAP{r}(idx,:));
                T.(cls{c}).run(end+1) = r;
            end
        end
        fprintf('%s: data size %s [samples x channels x trials]\n', cls{c}, mat2str(size(T.(cls{c}).mu)));
    end

    % Mu power per trial and channel: mean of squared samples over the trial
    for c = 1:2
        for f = {'mu','car','lap'}
            X = T.(cls{c}).(f{1});
            T.(cls{c}).(['P_' f{1}]) = squeeze(mean(X.^2, 1));   % [channels x trials]
        end
        fprintf('%s mu power size (no filter): %s\n', cls{c}, mat2str(size(T.(cls{c}).P_mu)));
    end

    % Last 0.5 s of each trial for the grand-average topoplot
    nL = round(lastWin_s*fs);
    for c = 1:2
        for f = {'mu','car','lap'}
            X = T.(cls{c}).(f{1})(end-nL+1:end, :, :);
            T.(cls{c}).(['GA_' f{1}]) = mean(squeeze(mean(X.^2, 1)), 2);   % [channels x 1]
        end
    end
    fprintf('Grand-average vectors ready for topoplot (32 x 1 each).\n');
end