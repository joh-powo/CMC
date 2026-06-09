% RUN_EXAMPLES  Run lightweight CMC encode/decode checks on the example set.
%
% This script auto-discovers every input listed in examples/manifest.csv and
% runs each one through the full CMC pipeline:
%
%     text -> Wubi -> Brotli + constrained DNA -> decode -> text
%
% For every example it reports basic sequence-quality metrics (DNA length,
% GC content, maximum homopolymer run, information density) and a round-trip
% PASS/FAIL check. The examples are intentionally small; they are meant for
% usability testing and demonstration, not for reproducing the manuscript
% benchmark numbers (which are based on data_raw/ and data_normalized/).
%
% Run from the project root:
%     run('examples/run_examples.m')

clc;

projectDir  = fileparts(fileparts(mfilename('fullpath')));
examplesDir = fullfile(projectDir, 'examples');

% --- Load Wubi-86 dictionary + Traditional extension table ---
opts = detectImportOptions(fullfile(projectDir, 'character.csv'));
opts.PreserveVariableNames = true;
wubiTable = readtable(fullfile(projectDir, 'character.csv'), opts);

optsTrad = detectImportOptions(fullfile(projectDir, 'trad_mapping.csv'));
optsTrad.PreserveVariableNames = true;
tradTable = readtable(fullfile(projectDir, 'trad_mapping.csv'), optsTrad);

mergedTable = [wubiTable; tradTable];
fprintf('Loaded Wubi table: %d rows; Traditional extension: %d rows.\n\n', ...
    height(wubiTable), height(tradTable));

% --- Discover example inputs from the manifest ---
manifestPath = fullfile(examplesDir, 'manifest.csv');
mopts = detectImportOptions(manifestPath);
mopts.PreserveVariableNames = true;
manifest = readtable(manifestPath, mopts);

exampleFiles = manifest.filename;
if ~iscell(exampleFiles)
    exampleFiles = cellstr(exampleFiles);
end
nFiles = numel(exampleFiles);

results   = cell(nFiles, 1);
passCount = 0;

for k = 1:nFiles
    fileName = strtrim(exampleFiles{k});
    filePath = fullfile(examplesDir, fileName);
    if exist(filePath, 'file') ~= 2
        fprintf('[%d/%d] %s  -- file not found, skipped\n', k, nFiles, fileName);
        continue;
    end

    text = strtrim(fileread(filePath));

    [wubi, punctMap, collisionMap] = wubi_encoder(text, mergedTable);
    [dna, meta]   = cmc_brotli_encoder(wubi);
    decodedWubi   = cmc_brotli_decoder(dna, meta);
    decodedText   = wubi_decoder(decodedWubi, mergedTable, punctMap, collisionMap);

    % --- Metrics ---
    nIn      = length(text);
    nWubi    = length(wubi);
    nDNA     = length(dna);
    gc       = 100 * sum(dna == 'G' | dna == 'C') / max(nDNA, 1);
    maxHomo  = maxHomopolymer(dna);
    infoBits = informationBits(text);
    density  = infoBits / max(nDNA, 1);
    ok       = strcmp(text, decodedText);
    if ok, passCount = passCount + 1; end

    results{k} = struct('file', fileName, 'in', nIn, 'wubi', nWubi, ...
        'dna', nDNA, 'gc', gc, 'homo', maxHomo, 'density', density, 'ok', ok);

    fprintf('\n[%d/%d] %s\n', k, nFiles, fileName);
    fprintf(['  input=%d chars | wubi=%d | dna=%d nt | GC=%.1f%% | ' ...
             'maxHomo=%d | density=%.2f bits/nt | %s\n'], ...
        nIn, nWubi, nDNA, gc, maxHomo, density, passFail(ok));
end

% --- Summary table ---
fprintf('\n========================= SUMMARY =========================\n');
fprintf('%-34s %5s %6s %6s %5s %8s  %s\n', ...
    'file', 'in', 'dna', 'GC%', 'homo', 'bits/nt', 'rt');
for k = 1:nFiles
    r = results{k};
    if isempty(r), continue; end
    fprintf('%-34s %5d %6d %6.1f %5d %8.2f  %s\n', ...
        r.file, r.in, r.dna, r.gc, r.homo, r.density, passFail(r.ok));
end
fprintf('-----------------------------------------------------------\n');
fprintf('Round-trip: %d / %d PASS\n', passCount, nFiles);

% ----------------------------- helpers -----------------------------
function label = passFail(ok)
    if ok
        label = 'PASS';
    else
        label = 'FAIL';
    end
end

function m = maxHomopolymer(dna)
    m = 0; run = 0; prev = char(0);
    for i = 1:length(dna)
        if dna(i) == prev
            run = run + 1;
        else
            run = 1; prev = dna(i);
        end
        if run > m, m = run; end
    end
end

function bits = informationBits(text)
    % Erlich-style content estimate: non-ASCII (CJK) = 16 bits, ASCII = 8 bits.
    bits = 0;
    for i = 1:length(text)
        if double(text(i)) > 127
            bits = bits + 16;
        else
            bits = bits + 8;
        end
    end
end
