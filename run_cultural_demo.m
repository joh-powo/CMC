clc;
clear all;

% 加载 OpenCC
opencc = py.importlib.import_module('opencc');
converter = opencc.OpenCC(pyargs('conversion', 't2s'));

% 加载码表
opts = detectImportOptions('character.csv');
opts.PreserveVariableNames = true;
wubiTable = readtable('character.csv', opts);
fprintf('加载码表: %d行\n', height(wubiTable));

% 繁体输入
text_trad = '國破山河在，城春草木深。龍飛鳳舞。';
fprintf('繁体输入: %s\n', text_trad);

% 繁转简
text_simp = char(converter.convert(text_trad));
fprintf('简体转换: %s\n', text_simp);

% 编码
wubi_code = wubi_encoder(text_simp, wubiTable);
fprintf('五笔码: %s (%d)\n', wubi_code, length(wubi_code));

[dna_seq, meta] = CMC_brotli_encoder(wubi_code);
fprintf('DNA: %d碱基, 密度: %.2f碱基/字\n', length(dna_seq), length(dna_seq)/length(text_trad));

% 解码
decoded_wubi = CMC_brotli_decoder(dna_seq, meta);
decoded_simp = wubi_decoder(decoded_wubi, wubiTable);
fprintf('解码(简): %s\n', decoded_simp);

% 简转繁还原
converter_s2t = opencc.OpenCC(pyargs('conversion', 's2t'));
decoded_trad = char(converter_s2t.convert(decoded_simp));
fprintf('还原(繁): %s\n', decoded_trad);

fprintf('匹配: %s\n', iif(strcmp(wubi_code, decoded_wubi), '✓', '✗'));

function r = iif(c, a, b)
    if c, r = a; else, r = b; end
end
