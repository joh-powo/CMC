clc;
clear all;

%% === 配置测试文件 ===
base_dir = fullfile(fileparts(mfilename('fullpath')), 'data_raw');
test_files = {
    fullfile(base_dir, '背影.docx'),         '简体存储';
    fullfile(base_dir, '繁体版出師表.docx'),  '繁体存简体(B方案)';
};

opencc = py.importlib.import_module('opencc');
t2s = opencc.OpenCC(pyargs('conversion', 't2s'));
opts = detectImportOptions('character.csv');
opts.PreserveVariableNames = true;
wubiTable = readtable('character.csv', opts);

results = struct([]);

for f = 1:size(test_files, 1)
    docx_file = test_files{f, 1};
    label = test_files{f, 2};
    fprintf('\n');
    fprintf('================================================================\n');
    fprintf('  测试 %d/%d: %s\n', f, size(test_files,1), label);
    fprintf('================================================================\n\n');

    %% === 读取文本 (空段落过滤) ===
    docx = py.importlib.import_module('docx');
    doc = docx.Document(docx_file);
    paragraphs = doc.paragraphs;

    text = '';
    for i = 0:int32(py.len(paragraphs))-1
        p = paragraphs{i+1};
        para_text = char(p.text);
        if ~isempty(strtrim(para_text))
            text = [text, para_text];
        end
    end
    fprintf('=== 输入 ===\n');
    fprintf('文件: %s\n', docx_file);
    fprintf('总字数: %d\n', length(text));
    fprintf('预览: %s...\n\n', text(1:min(50,end)));

    %% === 统一繁转简 (简体文本经 t2s 不变) ===
    text_simp = char(t2s.convert(text));
    fprintf('=== 简体化 ===\n');
    fprintf('预览: %s...\n\n', text_simp(1:min(50,end)));

    %% === 五笔编码 ===
    [wubi, punctMap, collisionMap] = wubi_encoder(text_simp, wubiTable);
    fprintf('=== 五笔编码 ===\n');
    fprintf('五笔码长度: %d\n', length(wubi));
    fprintf('标点数量: %d\n', length(punctMap.positions));
    fprintf('碰撞消歧: %d\n', length(collisionMap.positions));
    fprintf('预览: %s...\n\n', wubi(1:min(80,end)));

    %% === DNA 编码 ===
    [dna, meta] = cmc_brotli_encoder(wubi);
    meta.punctMap = punctMap;
    meta.collisionMap = collisionMap;
    fprintf('=== DNA编码 ===\n');
    fprintf('DNA长度: %d碱基\n', length(dna));
    fprintf('预览: %s...\n\n', dna(1:min(80,end)));

    %% === 统计 ===
    fprintf('=== 统计 ===\n');
    fprintf('原文字数: %d\n', length(text));
    fprintf('五笔码长: %d\n', meta.originalLen);
    fprintf('压缩后: %d 字节 (%.1f%%)\n', meta.compressedLen, 100*meta.compressedLen/meta.originalLen);
    fprintf('密度: %.2f 碱基/字\n', length(dna)/length(text));

    % Erlich公式: 汉字 16 bits, 英文 8 bits
    info_bits = 0; n_chinese = 0; n_english = 0;
    for i = 1:length(text)
        if double(text(i)) > 127
            info_bits = info_bits + 16;
            n_chinese = n_chinese + 1;
        else
            info_bits = info_bits + 8;
            n_english = n_english + 1;
        end
    end
    erlich_density = info_bits / length(dna);
    fprintf('存储密度(Erlich): %.2f bits/nt\n', erlich_density);
    fprintf('  其中汉字: %d (×16 bits), 英文/ASCII: %d (×8 bits)\n', n_chinese, n_english);

    gc = sum(dna == 'G' | dna == 'C');
    gc_pct = 100 * gc / length(dna);
    fprintf('GC含量: %.1f%%\n\n', gc_pct);

    %% === 解码验证 ===
    decoded_wubi = cmc_brotli_decoder(dna, meta);
    decoded_simp = wubi_decoder(decoded_wubi, wubiTable, meta.punctMap, meta.collisionMap);

    % 验证: 解码简体 vs t2s(原文)
    ref_simp = char(t2s.convert(text));
    wubi_ok = strcmp(wubi, decoded_wubi);
    text_ok = strcmp(ref_simp, decoded_simp);

    fprintf('=== 解码验证 ===\n');
    fprintf('解码字数: %d\n', length(decoded_simp));
    fprintf('解码预览: %s...\n', decoded_simp(1:min(50,end)));
    fprintf('参考简体: %s...\n', ref_simp(1:min(50,end)));
    fprintf('五笔匹配: %s\n', iif(wubi_ok, 'OK', 'FAIL'));
    fprintf('简体回环: %s\n', iif(text_ok, 'OK ✓✓✓', 'FAIL'));

    if ~text_ok
        minLen = min(length(ref_simp), length(decoded_simp));
        diffs = find(ref_simp(1:minLen) ~= decoded_simp(1:minLen));
        if ~isempty(diffs)
            fprintf('差异数量: %d\n', length(diffs));
            for d = 1:min(10, length(diffs))
                pos = diffs(d);
                fprintf('  位置 %d: 参考[%s](U+%04X) vs 解码[%s](U+%04X)\n', ...
                    pos, ref_simp(pos), double(ref_simp(pos)), decoded_simp(pos), double(decoded_simp(pos)));
            end
        end
    end

    %% === 不良基序检测 ===
    fprintf('\n=== 不良基序检测 ===\n');
    bad_motifs = {'TGC', 'CGC', 'GTC', 'GTG', 'GAC', 'CAC', 'GCG'};
    total_bad = 0;
    total_triplets = length(dna) - 2;
    for k = 1:length(bad_motifs)
        motif = bad_motifs{k};
        count = length(strfind(dna, motif));
        total_bad = total_bad + count;
        fprintf('  %s 出现次数: %d\n', motif, count);
    end
    bad_ratio = total_bad / total_triplets * 100;
    fprintf('不良基序总数: %d / %d 三联体\n', total_bad, total_triplets);
    fprintf('不良基序占比: %.3f%% (目标 < 0.6%%)\n', bad_ratio);
    fprintf('%s 不良基序约束\n', iif(bad_ratio < 0.6, '✓', '✗'));

    %% === 同聚物检测 ===
    fprintf('\n=== 同聚物检测 ===\n');
    max_run = 1; current_run = 1;
    for i = 2:length(dna)
        if dna(i) == dna(i-1)
            current_run = current_run + 1;
        else
            max_run = max(max_run, current_run);
            current_run = 1;
        end
    end
    max_run = max(max_run, current_run);
    fprintf('最大同聚物长度: %d (目标 ≤2)\n', max_run);

    %% === 自由能 & 汉明距离 ===
    fprintf('\n=== 自由能 & 汉明距离 ===\n');
    frag_len = 120;
    num_frags = min(6, floor(length(dna) / frag_len));
    fragments = cell(1, num_frags);
    for i = 1:num_frags
        start_idx = (i-1) * frag_len + 1;
        end_idx = start_idx + frag_len - 1;
        fragments{i} = dna(start_idx:end_idx);
    end
    fprintf('提取 %d 个 %dnt 片段用于检测\n\n', num_frags, frag_len);

    free_energies = zeros(1, num_frags);
    try
        RNA_mod = py.importlib.import_module('RNA');
        fprintf('使用 ViennaRNA 计算自由能...\n');
        for i = 1:num_frags
            seq_py = py.str(fragments{i});
            result = RNA_mod.fold(seq_py);
            mfe_val = double(result{2});
            free_energies(i) = mfe_val;
            fprintf('  片段 %d: ΔG = %.1f kcal/mol\n', i, mfe_val);
        end
    catch
        fprintf('ViennaRNA 不可用，使用简化近邻模型估算自由能...\n');
        nn_params = containers.Map();
        nn_params('AA') = -1.0; nn_params('TT') = -1.0;
        nn_params('AT') = -0.88; nn_params('TA') = -0.58;
        nn_params('CA') = -1.45; nn_params('TG') = -1.45;
        nn_params('GT') = -1.44; nn_params('AC') = -1.44;
        nn_params('CT') = -1.28; nn_params('AG') = -1.28;
        nn_params('GA') = -1.30; nn_params('TC') = -1.30;
        nn_params('CG') = -2.17; nn_params('GC') = -2.24;
        nn_params('GG') = -1.84; nn_params('CC') = -1.84;
        for i = 1:num_frags
            seq = fragments{i};
            dg = 0;
            for j = 1:length(seq)-1
                dinuc = seq(j:j+1);
                if isKey(nn_params, dinuc)
                    dg = dg + nn_params(dinuc);
                else
                    dg = dg - 1.0;
                end
            end
            free_energies(i) = dg;
            fprintf('  片段 %d: ΔG ≈ %.1f kcal/mol (近似估算)\n', i, dg);
        end
    end

    avg_fe = mean(free_energies);
    min_fe = min(free_energies);
    fprintf('\n平均自由能: %.1f kcal/mol\n', avg_fe);
    fprintf('最小自由能: %.1f kcal/mol (目标 > -30)\n', min_fe);
    fprintf('%s 自由能约束\n', iif(min_fe > -30, '✓', '✗'));

    fprintf('\n--- 汉明距离 ---\n');
    min_hamming = Inf;
    for i = 1:num_frags
        for j = i+1:num_frags
            hd = sum(fragments{i} ~= fragments{j});
            fprintf('  片段 %d vs %d: 汉明距离 = %d\n', i, j, hd);
            if hd < min_hamming, min_hamming = hd; end
        end
    end
    fprintf('\n最小汉明距离: %d (目标 > 60)\n', min_hamming);
    fprintf('%s 汉明距离约束\n', iif(min_hamming > 60, '✓', '✗'));

    %% === 保存结果 ===
    r.label = label;
    r.charCount = length(text);
    r.dnaLen = length(dna);
    r.density = erlich_density;
    r.gcPct = gc_pct;
    r.maxHomo = max_run;
    r.badRatio = bad_ratio;
    r.avgMFE = avg_fe;
    r.minHamming = min_hamming;
    r.textOK = text_ok;
    results = [results, r];
end

%% ========== 综合对比报告 ==========
fprintf('\n\n');
fprintf('╔══════════════════════════════════════════════════════════════════╗\n');
fprintf('║                CMC 编码质量综合对比报告                        ║\n');
fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
fprintf('║ %-20s', '指标');
for f = 1:length(results)
    fprintf(' │ %-16s', results(f).label);
end
fprintf(' ║\n');
fprintf('╠══════════════════════════════════════════════════════════════════╣\n');

fprintf('║ %-20s', '原文字数');
for f = 1:length(results), fprintf(' │ %12d 字 ', results(f).charCount); end
fprintf(' ║\n');

fprintf('║ %-20s', 'DNA碱基数');
for f = 1:length(results), fprintf(' │ %11d nt  ', results(f).dnaLen); end
fprintf(' ║\n');

fprintf('║ %-20s', '存储密度(Erlich)');
for f = 1:length(results), fprintf(' │ %10.2f bits/nt', results(f).density); end
fprintf(' ║\n');

fprintf('║ %-20s', 'GC含量');
for f = 1:length(results), fprintf(' │ %14.1f%%  ', results(f).gcPct); end
fprintf(' ║\n');

fprintf('║ %-20s', '最大同聚物');
for f = 1:length(results), fprintf(' │ %16d', results(f).maxHomo); end
fprintf(' ║\n');

fprintf('║ %-20s', '不良基序占比');
for f = 1:length(results), fprintf(' │ %13.3f%%  ', results(f).badRatio); end
fprintf(' ║\n');

fprintf('║ %-20s', '平均自由能');
for f = 1:length(results), fprintf(' │ %9.1f kcal/mol', results(f).avgMFE); end
fprintf(' ║\n');

fprintf('║ %-20s', '最小汉明距离');
for f = 1:length(results), fprintf(' │ %16d', results(f).minHamming); end
fprintf(' ║\n');

fprintf('║ %-20s', '全文匹配');
for f = 1:length(results), fprintf(' │ %16s', iif(results(f).textOK, 'OK ✓', 'FAIL')); end
fprintf(' ║\n');

fprintf('╚══════════════════════════════════════════════════════════════════╝\n');

function r = iif(c, a, b)
    if c, r = a; else, r = b; end
end
