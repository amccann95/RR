step = {'baseline', 'end ablation'};
coef = 0.1;

% [Ordbl,Ordend] = deal(zeros(41,1));

for m = 1:length(step)
    Valbl = zeros(41,5);
    St = zeros(100,5); 
    outcome = {};
    j = 1;
    liste = dir(strcat('data/',step{m},'/','*.mat'));
    pNb = 6;
    
    for k = pNb%1:41
        disp(cat(2,num2str(k),'  ',liste(k).name));
        filename = liste(k).name; load(strcat(step{m},'/',filename));
        outcome{k} = filename(13:end);
    
        RR=ecg_WS.RR;
        A = immersion(RR,2,1);
        [REC,DET,ENT,DIV,trend,tranche,RP] = recurrence_plot(A,coef,0);
        if (k == pNb)
           plot_rp_ecg(ecg_WS,REC,DET,RP,step{m}); 
        end
        Valbl(k,:) = [REC,DET,ENT,DIV,trend];
        St(1,:) = [REC,DET,ENT,DIV,trend];
        for p=2:100,
            sur = surrogate(RR,0);
            A = immersion(sur,2,1);
            [REC,DET,ENT,DIV,trend,tranche] = recurrence_plot(A,0.15,0);
            St(p,:) = [REC,DET,ENT,DIV,trend];
        end
        for p=1:5,
            [b,Ind]=sort(St(:,p),'descend');
            Ordbl(k,p) = find(Ind==1);
        end
         p = 2;
        [b,Ind] = sort(St(:,p),'descend');
        if (m == 1) 
            Ordbl(k) = find(Ind == 1);
        elseif (m == 2)
            Ordend(k) = find(Ind == 1);
        end
    end
% outcome = categorical(outcome);
% outcome = outcome';
% index = find(outcome == 'NLT_NoRec.mat');
% outcome(index) = [];
% Valbl(index,:) = [];
% Rec = Valbl(:,1);
% Det = Valbl(:,2);
% Ent = Valbl(:,3);
% Div = Valbl(:,4);
% trend = Valbl(:,5);
% 
% % compare LT vs NLT
% LTIndex = find((ismember(outcome,{'LT_NoRec.mat'})) == 1);
% LTIndex = [LTIndex; find((ismember(outcome,{'LT_Rec.mat'})) == 1)];
% NTIndex = find((ismember(outcome,{'NLT_Rec.mat'})) == 1);
% 
% recLTData(m,:) = Rec(LTIndex);
% recNTData(m,:) = Rec(NTIndex);
% detLTData(m,:) = Det(LTIndex);
% detNTData(m,:) = Det(NTIndex);

end

% Ordbl(33) = []; 
% Ordend(33) = [];

% pRecLT = permutationTest(recLTData(1,:),recLTData(2,:),1000,'sidedness','smaller');
% pRecNT = permutationTest(recNTData(1,:),recNTData(2,:),1000,'sidedness','smaller');
% 
% pDetLT = permutationTest(detLTData(1,:),detLTData(2,:),1000,'sidedness','smaller');
% pDetNT = permutationTest(detNTData(1,:),detNTData(2,:),1000,'sidedness','smaller');