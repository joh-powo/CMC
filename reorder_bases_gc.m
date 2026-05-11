function reordered = reorder_bases_gc(allowed, gc_count, at_count)
% REORDER_BASES_GC  Dynamically reorder allowed bases to balance GC content.
%   When running GC% > 50%, put AT bases first to reduce GC.
%   When running GC% < 50%, put GC bases first to increase GC.
%   When equal, use natural order.

    total = gc_count + at_count;
    if total == 0
        % No data yet, use natural order
        reordered = allowed;
        return;
    end
    
    gc_ratio = gc_count / total;
    
    is_gc = (allowed == 'G' | allowed == 'C');
    is_at = (allowed == 'A' | allowed == 'T');
    
    gc_bases = allowed(is_gc);
    at_bases = allowed(is_at);
    
    if gc_ratio > 0.5
        % GC too high, prefer AT first
        reordered = [at_bases, gc_bases];
    else
        % GC too low or equal, prefer GC first
        reordered = [gc_bases, at_bases];
    end
end
