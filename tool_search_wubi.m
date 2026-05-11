clc;
opts = detectImportOptions('character.csv');
opts.PreserveVariableNames = true;
T = readtable('character.csv', opts);

fprintf('Searching for "我"...\n');
for r = 1:height(T)
    for c = 2:min(5, width(T))
        v = T{r, c};
        if iscell(v), v = v{1}; end
        if ~isempty(v) && contains(v, '我')
            code = T{r, 1};
            if iscell(code), code = code{1}; end
            fprintf('Row %d: Code=%s, Char=%s\n', r, code, v);
        end
    end
end

fprintf('\nSearching for "actn"...\n');
for r = 1:height(T)
    code = T{r, 1};
    if iscell(code), code = code{1}; end
    if strcmp(code, 'actn')
        fprintf('Row %d: Code=actn, Chars=', r);
        for c = 2:min(5, width(T))
            v = T{r, c};
            if iscell(v), v = v{1}; end
            fprintf('%s ', v);
        end
        fprintf('\n');
    end
end