clc; clear all;

opts = detectImportOptions('character.csv');
opts.PreserveVariableNames = true;
wubiTable = readtable('character.csv', opts);

% 超长文本测试 (重复5遍)
base_text = ['我与父亲不相见已二年余了我最不能忘记的是他的背影那年冬天' ...
    '祖母死了父亲的差使也交卸了正是祸不单行的日子我从北京到徐州' ...
    '打算跟着父亲奔丧回家到徐州见着父亲看见满院狼藉的东西又想起祖母' ...
    '不禁簌簌地流下眼泪父亲说事已如此不必难过好在天无绝人之路'];

long_text = '';
for i = 1:5
    long_text = [long_text, base_text];
end

fprintf('原文: %d 字\n', length(long_text));

wubi_code = wubi_encoder(long_text, wubiTable);
fprintf('五笔码: %d\n', length(wubi_code));

[dna, meta] = CMC_brotli_encoder(wubi_code);
decoded = CMC_brotli_decoder(dna, meta);

fprintf('\n=== 结果 ===\n');
fprintf('DNA: %d 碱基\n', length(dna));
fprintf('密度: %.2f 碱基/字\n', length(dna)/length(long_text));

info_bits = length(long_text) * 13;
bits_per_nt = info_bits / length(dna);
fprintf('理论密度: %.2f bits/nt\n', bits_per_nt);
fprintf('匹配: %s\n', iif(strcmp(wubi_code, decoded), '✓', '✗'));

function r = iif(c, a, b)
    if c, r = a; else, r = b; end
end
