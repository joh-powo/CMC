function text = wubi_decoder(wubi_str, wubiTable, punctMap, collisionMap)
    % Build code -> ordered char list (same traversal order as encoder)
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

    % Parse collision map into position-indexed lookup
    hasCollisionMap = nargin >= 4 && ~isempty(collisionMap) && ...
                      ~isempty(collisionMap.positions);
    collisionLookup = containers.Map('KeyType', 'int32', 'ValueType', 'int32');
    if hasCollisionMap
        for k = 1:length(collisionMap.positions)
            collisionLookup(int32(collisionMap.positions(k))) = ...
                int32(collisionMap.indices(k));
        end
    end

    % Decode wubi segments to characters
    segments = strsplit(wubi_str, '6');
    chars = {};
    for k = 1:length(segments)
        seg = segments{k};
        if isKey(codeToChars, seg)
            candidates = codeToChars(seg);
            if length(candidates) > 1 && isKey(collisionLookup, int32(k))
                % Use collision map to pick the correct character
                idx = collisionLookup(int32(k));
                if idx >= 1 && idx <= length(candidates)
                    chars{end+1} = candidates(idx);
                else
                    chars{end+1} = candidates(1);
                end
            else
                % No collision or no map entry: use first candidate
                chars{end+1} = candidates(1);
            end
        end
    end
    
    % Re-insert punctuation
    if nargin >= 3 && ~isempty(punctMap) && ~isempty(punctMap.positions)
        text = '';
        punctIdx = 1;
        for i = 1:length(chars)
            while punctIdx <= length(punctMap.positions) && punctMap.positions(punctIdx) < i
                text = [text, punctMap.chars{punctIdx}];
                punctIdx = punctIdx + 1;
            end
            text = [text, chars{i}];
        end
        while punctIdx <= length(punctMap.positions)
            text = [text, punctMap.chars{punctIdx}];
            punctIdx = punctIdx + 1;
        end
    else
        text = strjoin(chars, '');
    end
    
    fprintf('Decoded %d chars.\n', length(text));
end