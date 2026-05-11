% export_dna_for_plots.m - Export DNA sequences and stage data for plotting
clc; clear all;

%% === Traditional Chinese (Chu Shi Biao) ===
docx_file = 'E:\研\李孜淇\Projects\CMC附件\繁体版出師表.docx';
docx = py.importlib.import_module('docx');
opencc = py.importlib.import_module('opencc');
t2s = opencc.OpenCC(pyargs('conversion', 't2s'));
opts = detectImportOptions('character.csv');
opts.PreserveVariableNames = true;
wubiTable = readtable('character.csv', opts);

doc = docx.Document(docx_file);
paragraphs = doc.paragraphs;
text = '';
for i = 0:int32(py.len(paragraphs))-1
    p = paragraphs{i+1};
    text = [text, char(p.text)];
end
text_simp = char(t2s.convert(text));
[wubi, ~, ~] = wubi_encoder(text_simp, wubiTable);
[dna, meta] = CMC_brotli_encoder(wubi);

fid = fopen('plots/dna_trad.txt', 'w');
fprintf(fid, '%s', dna);
fclose(fid);
fprintf('Trad: %d chars, %d wubi, %d compressed, %d nt\n', ...
    length(text), length(wubi), meta.compressedLen, length(dna));

%% === Simplified Chinese (Beiying) ===
docx_file2 = 'E:\研\李孜淇\Projects\CMC附件\背影test.docx';
doc2 = docx.Document(docx_file2);
paragraphs2 = doc2.paragraphs;
text2 = '';
for i = 0:int32(py.len(paragraphs2))-1
    p = paragraphs2{i+1};
    para_text = char(p.text);
    if ~isempty(strtrim(para_text))
        text2 = [text2, para_text];
    end
end
[wubi2, ~, ~] = wubi_encoder(text2, wubiTable);
[dna2, meta2] = CMC_brotli_encoder(wubi2);

fid = fopen('plots/dna_simp.txt', 'w');
fprintf(fid, '%s', dna2);
fclose(fid);
fprintf('Simp: %d chars, %d wubi, %d compressed, %d nt\n', ...
    length(text2), length(wubi2), meta2.compressedLen, length(dna2));

fprintf('Export done.\n');

