function [wubiCode, punctMap, collisionMap] = wubi_encoder(inputStr, wubiTable)
    % Build char -> code mapping
    charToCode = containers.Map;
    for row = 1:height(wubiTable)
        code = wubiTable{row, 1};
        if iscell(code), code = code{1}; end
        for col = 2:min(5, width(wubiTable))
            charVal = wubiTable{row, col};
            if iscell(charVal), charVal = charVal{1}; end
            if ~isempty(charVal) && ischar(charVal)
                if length(charVal) == 1
                    charToCode(charVal) = code;
                end
            end
        end
    end
    fprintf('Built lookup map with %d characters.\n', length(keys(charToCode)));

    % Build code -> ordered char list (for collision detection)
    % Traversal order MUST match wubi_decoder exactly
    codeToChars = containers.Map;
    for row = 1:height(wubiTable)
        code = wubiTable{row, 1};
        if iscell(code), code = code{1}; end
        for col = 2:min(5, width(wubiTable))
            charVal = wubiTable{row, col};
            if iscell(charVal), charVal = charVal{1}; end
            if ~isempty(charVal) && ischar(charVal) && length(charVal) == 1
                if isKey(codeToChars, code)
                    existing = codeToChars(code);
                    if ~contains(existing, charVal)
                        codeToChars(code) = [existing, charVal];
                    end
                else
                    codeToChars(code) = charVal;
                end
            end
        end
    end

    wubiCode = '';
    punctMap = struct('positions', [], 'chars', {{}});
    collisionMap = struct('positions', [], 'indices', []);
    charIdx = 0;
    
    for i = 1:length(inputStr)
        ch = inputStr(i);
        if isKey(charToCode, ch)
            code = charToCode(ch);
            if ~isempty(wubiCode)
                wubiCode = [wubiCode, '6'];
            end
            wubiCode = [wubiCode, code];
            charIdx = charIdx + 1;

            % Collision detection: if this code maps to >1 char, record index
            if isKey(codeToChars, code)
                candidates = codeToChars(code);
                if length(candidates) > 1
                    idx = find(candidates == ch, 1);
                    if isempty(idx), idx = 1; end
                    collisionMap.positions(end+1) = charIdx;
                    collisionMap.indices(end+1) = idx;
                end
            end
        else
            punctMap.positions(end+1) = charIdx;
            punctMap.chars{end+1} = ch;
        end
    end
    fprintf('Skipped %d punctuation chars.\n', length(punctMap.positions));
    fprintf('Recorded %d collision disambiguations.\n', length(collisionMap.positions));
end