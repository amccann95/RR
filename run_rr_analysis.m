% run rr_analysis.m
taus = [1];
ms = [2];
coefs = [0.1, 0.15];
steps = {'baseline', 'end ablation'};

for step = steps
    for tau = taus
        for m = ms
            for coef = coefs
                rr_analysis(step{1},tau,m,coef,pwd);
            end
        end
    end
end