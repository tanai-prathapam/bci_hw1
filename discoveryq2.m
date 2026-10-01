rng(1);
dp = fullfile('ErrP_data_scripts','ErrP_data_scripts');

%% 1. Data file
whos('-file', fullfile(dp,'ErrP_data_HW1.mat'))
D = load(fullfile(dp,'ErrP_data_HW1.mat'));
assert(isfield(D,'trainingEpochs'), 'Expected top-level variable trainingEpochs');
E = D.trainingEpochs;
disp(fieldnames(E)')

X  = E.rotation_data;
sz = size(X), cls = class(X)
assert(ndims(X)==3 && sz(1)==1024 && sz(2)==32, 'Expected [time x ch x trials] = [1024 32 N]');
fs = sz(1)/2                                  % expect 512
t  = -1 + (0:sz(1)-1)'/fs;                    % t(513) = 0
fprintf('t(1)=%.4f  t(513)=%.4f  t(end)=%.4f\n', t(1), t(513), t(end));

%% 2. Labels, magnitudes, counts (replaces crosstab)
labs = unique(E.label(:)); mags = unique(E.magnitude(:));
[~,li] = ismember(E.label(:), labs); [~,mi] = ismember(E.magnitude(:), mags);
disp(array2table(accumarray([li mi],1), ...
    'RowNames', strcat("label",string(labs)), 'VariableNames', strcat("mag",string(mags))))
fileIDs = unique(E.fileID(:))', sessionIDs = unique(E.sessionID(:))'

%% 3. NaN and scale checks (replaces prctile)
hasNaN = any(isnan(X(:)))
v = sort(abs(X(:)));
scale_50_99_99p9 = v(max(1,round([0.50 0.99 0.999]*numel(v))))'

%% 4. Channel file and automatic Cz / FCz lookup
C  = load(fullfile(dp,'ErrP_channels.mat'));
P  = C.params;
CL = P.chanlocs;
disp(fieldnames(CL)')
disp(CL(1))

lf = intersect({'labels','label','name'}, fieldnames(CL), 'stable');
assert(~isempty(lf), 'No label field in chanlocs; check the printout above');
chanNames = strtrim(arrayfun(@(c) char(c.(lf{1})), CL, 'UniformOutput', false));
disp(chanNames)

czIdx  = find(strcmpi(chanNames,'Cz'));
fczIdx = find(strcmpi(chanNames,'FCz'));
assert(isscalar(czIdx), 'Cz not found exactly once; read chanNames above');
fprintf('Cz = %d, FCz = %s\n', czIdx, mat2str(fczIdx));

t  = P.epochTime(:);
fs = P.fsamp;
assert(abs(t(P.epochOnset)) < 1e-6, 'epochOnset is not t = 0');

%% 5. Topoplot header check (first 6 lines)
fid = fopen(fullfile(dp,'topoplot.m'));
if fid > 0
    for k = 1:6, disp(fgetl(fid)); end, fclose(fid);
else
    fprintf('topoplot.m not found in %s (check the project root)\n', dp);
end