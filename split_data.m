% script to split patient files

step = 'end ablation';
rd_dir = strcat('./data/', step, '/');
wr_dir = strcat('./data_split/',step,'/');
pts = [1:41];
n_pts = length(pts);
n_segs_file = 5; % seconds

for pt = 1:n_pts
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
    n_samples = length(ecg);
    fs = data_struct.ecg_WS.Samplingfrequency;

    n_samps_file = floor(n_samples / n_segs_file);

    start_idx = 1;
    end_idx = n_samps_file;
    for file_nb = 1:n_segs_file
        split_ecg = ecg(start_idx:end_idx);
        rr_ind_split = rr_ind(rr_ind >= start_idx & rr_ind <= end_idx);
        ecg_WS.R_inds = rr_ind_split;
        ecg_WS.RR = diff(rr_ind_split) / 2;  % convert from samples to ms
        ecg_WS.data = split_ecg;

        fname = strcat(pt_char,'-',num2str(file_nb));
        save(strcat(wr_dir,fname),'ecg_WS');

        start_idx = start_idx + n_samps_file;
        end_idx = end_idx + n_samps_file;
    end
end