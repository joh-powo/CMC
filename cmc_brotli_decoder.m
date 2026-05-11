function decoded_str = cmc_brotli_decoder(dna_sequence, meta)
% CMC_BROTLI_DECODER  Decode constrained DNA sequence back to wubi string.
%   Reverses the encoding with limited ternary grouping.

    dna_len = length(dna_sequence);
    nbits = meta.totalBits;
    max_ternary_groups = meta.maxTernaryGroups;

    % --- Step 1: DNA to bitstream ---
    bitstream = zeros(1, nbits);
    gc_count = 0;
    at_count = 0;

    % Track encoder's bit consumption for position-aware writing
    bits_consumed = 0;
    ternary_groups_used = 0;

    % Ternary group state
    ternary_buf = [];
    ternary_write_pos = 0;

    for i = 1:dna_len
        if i < 3
            allowed = 'ACGT';
        else
            allowed = get_allowed_bases(dna_sequence(i-2:i-1));
        end

        allowed = reorder_bases_gc(allowed, gc_count, at_count);
        n = length(allowed);
        idx = find(allowed == dna_sequence(i)) - 1;

        if n >= 4
            bits_remaining = nbits - bits_consumed;
            if bits_remaining >= 2
                bitstream(bits_consumed + 1) = floor(idx/2);
                bitstream(bits_consumed + 2) = mod(idx, 2);
                bits_consumed = bits_consumed + 2;
            elseif bits_remaining == 1
                bitstream(bits_consumed + 1) = idx;
                bits_consumed = bits_consumed + 1;
            end
        elseif n == 3
            if ~isempty(ternary_buf)
                % Continue collecting ternary digits
                ternary_buf = [ternary_buf, idx];
                if length(ternary_buf) == 7
                    val = 0;
                    for t = 1:7
                        val = val * 3 + ternary_buf(t);
                    end
                    for b = 0:10
                        bitstream(ternary_write_pos + b + 1) = ...
                            bitand(floor(val / 2^(10-b)), 1);
                    end
                    ternary_buf = [];
                end
            else
                % Decision: ternary group or 1-bit?
                bits_remaining = nbits - bits_consumed;
                if ternary_groups_used < max_ternary_groups && ...
                   bits_remaining >= 30
                    % Encoder started a ternary group
                    ternary_write_pos = bits_consumed;
                    bits_consumed = bits_consumed + 11;
                    ternary_buf = [idx];
                    ternary_groups_used = ternary_groups_used + 1;
                elseif bits_remaining >= 1
                    % 1-bit fallback
                    bitstream(bits_consumed + 1) = idx;
                    bits_consumed = bits_consumed + 1;
                end
            end
        else
            bits_remaining = nbits - bits_consumed;
            if bits_remaining >= 1
                bitstream(bits_consumed + 1) = idx;
                bits_consumed = bits_consumed + 1;
            end
        end

        if dna_sequence(i) == 'G' || dna_sequence(i) == 'C'
            gc_count = gc_count + 1;
        else
            at_count = at_count + 1;
        end
    end

    % Handle remaining ternary digits
    if ~isempty(ternary_buf)
        nrem = length(ternary_buf);
        val = 0;
        for t = 1:nrem
            val = val * 3 + ternary_buf(t);
        end
        val = val * 3^(7 - nrem);
        for b = 0:10
            pos = ternary_write_pos + b + 1;
            if pos <= nbits
                bitstream(pos) = bitand(floor(val / 2^(10-b)), 1);
            end
        end
    end

    % --- Step 2: Bitstream to bytes ---
    nbytes = meta.compressedLen;
    bytes = zeros(1, nbytes);
    for i = 1:nbytes
        v = 0;
        for j = 0:7
            bidx = (i-1)*8 + j + 1;
            if bidx <= length(bitstream)
                v = bitor(bitshift(v, 1), bitstream(bidx));
            else
                v = bitshift(v, 1);
            end
        end
        bytes(i) = v;
    end

    % --- Step 3: Brotli decompress ---
    brotli = py.importlib.import_module('brotli');
    py_bytes = py.bytes(uint8(bytes));
    decompressed = brotli.decompress(py_bytes);
    decoded_str = char(uint8(py.array.array('B', decompressed)));

    fprintf('Decoded: %d bases -> %d symbols\n', dna_len, length(decoded_str));
end