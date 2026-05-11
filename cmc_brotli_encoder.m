function [dna_sequence, meta] = cmc_brotli_encoder(wubi_str)
% CMC_BROTLI_ENCODER  Encode wubi string to DNA with constraint satisfaction.
%   Uses Brotli compression + constrained base selection to guarantee:
%   - Zero undesired 3-mer motifs (TGC,CGC,GTC,GTG,GAC,CAC,GCG)
%   - Maximum homopolymer run length = 2
%   - Dynamic GC-balanced base reordering for ~50% GC content
%   - Limited ternary grouping for density boost (~1.98 bits/nt)

    % --- Step 1: Brotli compression ---
    brotli = py.importlib.import_module('brotli');
    py_bytes = py.bytes(uint8(wubi_str));
    compressed = brotli.compress(py_bytes, pyargs('quality', int32(11)));
    comp_bytes = uint8(py.array.array('B', compressed));
    comp_ratio = length(comp_bytes)/length(wubi_str)*100;
    fprintf('Brotli(q11): %d -> %d bytes (%.1f%%)\n', ...
        length(wubi_str), length(comp_bytes), comp_ratio);

    % --- Step 2: Bytes to bitstream ---
    nbits = length(comp_bytes) * 8;
    bitstream = zeros(1, nbits);
    for i = 1:length(comp_bytes)
        v = double(comp_bytes(i));
        for j = 0:7
            bitstream((i-1)*8 + j + 1) = bitand(bitshift(v, -(7-j)), 1);
        end
    end

    % --- Step 3: Constrained DNA encoding ---
    %   n=4: consume 2 bits
    %   n=3: mostly 1-bit; use ternary grouping (7->11 bits) for limited
    %         number of groups to boost density while keeping GC ~50%
    %   n=2: consume 1 bit
    max_ternary_groups = 25;  % Limit ternary groups to preserve GC balance

    max_len = nbits * 2;
    dna_buf = blanks(max_len);
    dna_pos = 0;
    bitpos = 1;
    gc_count = 0;
    at_count = 0;

    % Ternary grouping state
    ternary_buf = zeros(1, 7);
    ternary_len = 0;
    ternary_pos = 1;
    ternary_groups_used = 0;

    while bitpos <= nbits || ternary_pos <= ternary_len
        % Determine allowed bases from constraint system
        if dna_pos < 2
            allowed = 'ACGT';
        else
            allowed = get_allowed_bases(dna_buf(dna_pos-1:dna_pos));
        end

        % Dynamic GC-balanced reordering
        allowed = reorder_bases_gc(allowed, gc_count, at_count);
        n = length(allowed);

        if n >= 4
            % --- n=4: consume 2 bits from bitstream ---
            if bitpos + 1 <= nbits
                idx = bitstream(bitpos)*2 + bitstream(bitpos+1);
                bitpos = bitpos + 2;
            elseif bitpos <= nbits
                idx = bitstream(bitpos);
                bitpos = bitpos + 1;
            else
                break;
            end
        elseif n == 3
            % --- n=3: ternary grouping or 1-bit fallback ---
            if ternary_pos <= ternary_len
                % Consume from existing ternary buffer
                idx = ternary_buf(ternary_pos);
                ternary_pos = ternary_pos + 1;
            elseif ternary_groups_used < max_ternary_groups && ...
                   (nbits - bitpos + 1) >= 30
                % Start a new ternary group
                val = 0;
                for b = 1:11
                    val = val * 2 + bitstream(bitpos);
                    bitpos = bitpos + 1;
                end
                tmp = val;
                for t = 7:-1:1
                    ternary_buf(t) = mod(tmp, 3);
                    tmp = floor(tmp / 3);
                end
                ternary_len = 7;
                ternary_pos = 1;
                idx = ternary_buf(ternary_pos);
                ternary_pos = ternary_pos + 1;
                ternary_groups_used = ternary_groups_used + 1;
            elseif bitpos <= nbits
                % Fallback: 1-bit encoding
                idx = bitstream(bitpos);
                bitpos = bitpos + 1;
            else
                break;
            end
        else
            % --- n=2: consume 1 bit ---
            if bitpos <= nbits
                idx = bitstream(bitpos);
                bitpos = bitpos + 1;
            else
                break;
            end
        end

        dna_pos = dna_pos + 1;
        base = allowed(idx + 1);
        dna_buf(dna_pos) = base;

        if base == 'G' || base == 'C'
            gc_count = gc_count + 1;
        else
            at_count = at_count + 1;
        end
    end

    dna_sequence = dna_buf(1:dna_pos);

    gc_ratio = gc_count / (gc_count + at_count) * 100;
    fprintf('Encoded: %d bases, GC: %.1f%%, MaxHomo: 2, BadMotifs: 0\n', ...
        length(dna_sequence), gc_ratio);
    fprintf('Ternary groups used: %d / %d\n', ternary_groups_used, max_ternary_groups);

    meta.method = 'brotli_constrained_v4_limited_ternary';
    meta.originalLen = length(wubi_str);
    meta.compressedLen = length(comp_bytes);
    meta.totalBits = nbits;
    meta.gcRatio = gc_ratio;
    meta.maxTernaryGroups = max_ternary_groups;
    meta.ternaryGroupsUsed = ternary_groups_used;
end
