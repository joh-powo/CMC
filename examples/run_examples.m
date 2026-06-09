% RUN_EXAMPLES  Run lightweight CMC checks on the example input set.
% Execute from the project root:
%   run('examples/run_examples.m')

clc;

projectDir = fileparts(fileparts(mfilename('fullpath')));
examplesDir = fullfile(projectDir, 'examples');

exampleFiles = {
    'quick_test.txt'
    'simplified_modern.txt'
    'traditional_classical.txt'
    'punctuation_mixed.txt'
    'mixed_alnum_symbols.txt'
    'long_simplified_paragraph.txt'
    'traditional_variants.txt'
};

opts = detectImportOptions(fullfile(projectDir, 'character.csv'));
opts.PreserveVariableNames = true;
wubiTable = readtable(fullfile(projectDir, 'character.csv'), opts);

optsTrad = detectImportOptions(fullfile(projectDir, 'trad_mapping.csv'));
optsTrad.PreserveVariableNames = true;
tradTable = readtable(fullfile(projectDir, 'trad_mapping.csv'), optsTrad);
mergedTable = [wubiTable; tradTable];

fprintf('Loaded Wubi table: %d rows; Traditional extension: %d rows.\n', ...
    height(wubiTable), height(tradTable));

for k = 1:numel(exampleFiles)
    fileName = exampleFiles{k};
    filePath = fullfile(examplesDir, fileName);
    text = strtrim(fileread(filePath));

    fprintf('\n[%d/%d] %s\n', k, numel(exampleFiles), fileName);
    fprintf('Input length: %d chars\n', length(text));

    [wubi, punctMap, collisionMap] = wubi_encoder(text, mergedTable);
    [dna, meta] = cmc_brotli_encoder(wubi);
    decodedWubi = cmc_brotli_decoder(dna, meta);
    decodedText = wubi_decoder(decodedWubi, mergedTable, punctMap, collisionMap);

    gc = 100 * sum(dna == 'G' | dna == 'C') / length(dna);
    roundTripOK = strcmp(text, decodedText);

    fprintf('Wubi length: %d symbols\n', length(wubi));
    fprintf('DNA length: %d nt\n', length(dna));
    fprintf('GC content: %.1f%%\n', gc);
    fprintf('Round trip: %s\n', passFail(roundTripOK));
end

function label = passFail(ok)
    if ok
        label = 'PASS';
    else
        label = 'FAIL';
    end
end
