function rr_analysis(step, tau, m, coef, project_path)
% calculate time, frequency, and nonlinear indices for RR intervals. 

cd(project_path);

% parameters
pts = [1:32,34:41];
n_pts = length(pts);
rd_dir = strcat('./data/',step,'/');

% preallocate variables for indices
[rec, det, ent, div, sampen] = deal(zeros(n_pts,1));
[pnn20, pnn50, sdnn, rmssd] = deal(zeros(n_pts,1));

% calculate values of indices for each patient
idx = 1;
for pt = pts
    % calculate indices for all patient files
    if pt < 10
        pt_char = strcat('0',num2str(pt));
    else
        pt_char = num2str(pt);
    end

    % load in data
    data_struct = load(strcat(rd_dir,pt_char));
    rr_ind = data_struct.ecg_WS.R_inds;
    rr_int = data_struct.ecg_WS.RR;
    ecg = data_struct.ecg_WS.data;
    
    % get matrix of state-space vectors
    emb_vec = embedding(rr_int,tau,m);
    
    % calculate recurrence/nonlinear measures, of irregularity
    [rec(idx), det(idx), ent(idx), div(idx)] = recurrence_plot(emb_vec,coef,0);

    sampen(idx) = sample_entropy(rr_int,m,coef,'chebychev');

    % calculate time-domain measures of variability
    [pnn20(idx), pnn50(idx), sdnn(idx), rmssd(idx)] = time_var_measures(rr_int);
    
    idx = idx + 1;
end

rr_indices_struct.rec = rec;
rr_indices_struct.det = det;
rr_indices_struct.div = div;
rr_indices_struct.sampen = sampen;
rr_indices_struct.pnn20 = pnn20;
rr_indices_struct.pnn50 = pnn50;
rr_indices_struct.sdnn = sdnn;
rr_indices_struct.rmssd = rmssd;

% save rr_int indices
write_dir = strcat('./data/rr_indices/',step,'/');
if ~exist(write_dir,'dir')
    mkdir(write_dir);
end
write_filename = strcat('tau',num2str(tau),'_dim',num2str(m),'_coef',num2str(coef));
save(strcat(write_dir, write_filename,'.mat'),'rr_indices_struct');

end