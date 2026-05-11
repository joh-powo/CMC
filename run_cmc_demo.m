clc;
clear all;

% 1. 加载五笔表
fprintf('加载五笔表...\n');
opts = detectImportOptions('character.csv');
opts.PreserveVariableNames = true;
wubiTable = readtable('character.csv', opts);
fprintf('已加载 %d 行。\n', height(wubiTable));

% 2. 输入中文文本
original_text = '我与父亲不相见已二年余了';
fprintf('\n原文: %s\n', original_text);
fprintf('字数: %d\n', length(original_text));

% 3. 五笔编码 (Part A)
fprintf('\n[Part A] 五笔编码...\n');
wubi_code = wubi_encoder(original_text, wubiTable);
fprintf('五笔码: %s\n', wubi_code);
fprintf('码长: %d\n', length(wubi_code));

% 4. 正负码表编码 (Part B)
fprintf('\n[Part B] 正负码表编码...\n');
dna_seq = cmc_dual_table_encoder(wubi_code);

% 5. 结果展示
fprintf('\n========== 编码结果 ==========\n');
fprintf('DNA序列 (前60碱基): %s...\n', dna_seq(1:min(60, length(dna_seq))));
fprintf('总碱基数: %d\n', length(dna_seq));

% 6. 约束验证
fprintf('\n========== 约束验证 ==========\n');
max_run = 1;
current_run = 1;
for i = 2:length(dna_seq)
    if dna_seq(i) == dna_seq(i-1)
        current_run = current_run + 1;
    else
        max_run = max(max_run, current_run);
        current_run = 1;
    end
end
max_run = max(max_run, current_run);
fprintf('最大同聚物长度: %d (目标 ≤3)\n', max_run);

gc_count = sum(dna_seq == 'G') + sum(dna_seq == 'C');
gc_content = (gc_count / length(dna_seq)) * 100;
fprintf('GC含量: %.2f%% (目标 45-55%%)\n', gc_content);

% 7. 密度分析
fprintf('\n========== 密度分析 ==========\n');
density = length(dna_seq) / length(original_text);
fprintf('存储密度: %.2f 碱基/字\n', density);

% ========== 解码测试 ==========
fprintf('\n========== 解码测试 ==========\n');
decoded_wubi = cmc_dual_table_decoder(dna_seq);
fprintf('解码五笔: %s\n', decoded_wubi);

decoded_text = wubi_decoder(decoded_wubi, wubiTable);
fprintf('解码文本: %s\n', decoded_text);

if strcmp(wubi_code, decoded_wubi)
    fprintf('✓ 五笔码完全匹配\n');
else
    fprintf('✗ 五笔码不匹配\n');
end
