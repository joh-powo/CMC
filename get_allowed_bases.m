function allowed = get_allowed_bases(prev2)
% GET_ALLOWED_BASES  Return allowed next bases given 2-base context.
%   Enforces: (1) no undesired 3-mer motifs, (2) homopolymer run <= 2.
    persistent lookup;
    if isempty(lookup)
        lookup = containers.Map;
        bases = 'ACGT';
        bad_motifs = {'TGC','CGC','GTC','GTG','GAC','CAC','GCG'};
        for i = 1:4
            for j = 1:4
                ctx = [bases(i), bases(j)];
                forbidden = '';
                % Check each bad motif
                for m = 1:length(bad_motifs)
                    mot = bad_motifs{m};
                    if ctx(1)==mot(1) && ctx(2)==mot(2)
                        forbidden = [forbidden, mot(3)];
                    end
                end
                % Homopolymer: if last 2 are same, forbid that base
                if ctx(1) == ctx(2)
                    forbidden = [forbidden, ctx(1)];
                end
                forbidden = unique(forbidden);
                lookup(ctx) = bases(~ismember(bases, forbidden));
            end
        end
    end
    allowed = lookup(prev2);
end
