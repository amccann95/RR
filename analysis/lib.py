import scipy.io as spio
import numpy as np
import pandas as pd
from sklearn.metrics import auc
import matplotlib.pyplot as plt
from sklearn.model_selection import ShuffleSplit
from sklearn.feature_selection import f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve
from sklearn.metrics import RocCurveDisplay

def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict        

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict


def plot_mean_roc(fig, ax, tprs, aucs, mean_fpr, thresholds, var_name, title_string):
    
    ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(
        mean_fpr,
        mean_tpr,
        color="b",
        label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
        lw=2,
        alpha=0.8,
    )
    
    J = mean_tpr - mean_fpr
    idx = np.argmax(J)
    opt_th = thresholds[idx]
    sens = mean_tpr[idx]
    spec = 1 - mean_fpr[idx]
    
    plt.scatter(mean_fpr[idx], mean_tpr[idx], marker='o', color='black', label='Best')

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(
        mean_fpr,
        tprs_lower,
        tprs_upper,
        color="grey",
        alpha=0.2,
        label=r"$\pm$ 1 std. dev.",
    )

    ax.set(
        xlim=[-0.05, 1.05],
        ylim=[-0.05, 1.05],
        title=title_string,
    )
    ax.legend(loc="lower right")
    
    return plt, sens, spec, opt_th


def retrieve_df(tau, m, coef, feat_name, outcomes, outcome_names, pt_nbs):
    ROOT = "C:/Users/amcca/switchdrive/PhD/RR/"
    df = []
 
    # load in matlab data with appropriate params
    filename = 'tau' + str(tau) + '_dim' + str(m) + '_coef' + str(coef) + '.mat'
    mat_bl = loadmat(ROOT + 'data/rr_indices/baseline/' + filename)
    mat_end = loadmat(ROOT + 'data/rr_indices/end ablation/' + filename)

    # get each AF complexity index in column, each row corresponds to one patient only
    data_dict_bl = mat_bl['rr_indices_struct']
    data_dict_end = mat_end['rr_indices_struct']

    # append pandas data frame entry with all parameters
    for pt_nb in pt_nbs:
        df.append(
            {
                feat_name: data_dict_bl[feat_name][pt_nb],
                'tau': tau,
                'm': m,
                'coef': coef,
                'outcome': outcome_names[outcomes[pt_nb]],
                'step': 'BL'
            })

        df.append(
            {
                feat_name: data_dict_end[feat_name][pt_nb],
                'tau': tau,
                'm': m,
                'coef': coef,
                'outcome': outcome_names[outcomes[pt_nb]],
                'step': 'end_ABL'
            })
                    
     # convert dictionary to pandas dataframe for use with seaborn
    df = pd.DataFrame(df)
    
    return df


def get_p_vals(taus, ms, coefs, feat_name, y, step):
    ROOT = "C:/Users/amcca/switchdrive/PhD/RR/"
    df = []
    
    for tau in taus:
        for m in ms:
            for coef in coefs:
                # load in matlab data with appropriate params
                filename = 'tau' + str(tau) + '_dim' + str(m) + '_coef' + str(coef) + '.mat'
                if step == 'baseline':
                    mat = loadmat(ROOT + 'data/rr_indices/baseline/' + filename)
                elif step == 'end ablation':
                    mat = loadmat(ROOT + 'data/rr_indices/end ablation/' + filename)
                    
                data_dict = mat['rr_indices_struct']
                idx = list(data_dict).index(feat_name)
                data_matrix = np.stack(list(data_dict.values()))
                data_matrix = data_matrix.transpose()
                
                f = []
                p = []
                feat_names = data_dict.keys()
                p_norm = np.zeros((3, len(feat_names)))
                ss = ShuffleSplit(n_splits=3, test_size=0.2, random_state=0)
                
                for i, (train, test) in enumerate(ss.split(data_matrix)):
                    f_fold, p_fold = f_classif(data_matrix[train,idx].reshape(-1, 1), y[train])
                    f.append(p_fold)
                    p.append(p_fold)
                    
                f = np.vstack(f)
                p = np.vstack(p)
                
                p_val = np.mean(p, axis=0)
                
                df.append(
                    {
                        'tau': tau,
                        'm': m,
                        'coef': coef,
                        'p_val': p_val
                    })
                              
    df = pd.DataFrame(df)
    
    return df


def get_pred_stats(taus, ms, coefs, feat_name, y, step):
    ROOT = "C:/Users/amcca/switchdrive/PhD/RR/"
    df = []    
             
    for tau in taus:
        for m in ms:
            for coef in coefs:
                # load in matlab data with appropriate params
                filename = 'tau' + str(tau) + '_dim' + str(m) + '_coef' + str(coef) + '.mat'
                if step == 'baseline':
                    mat = loadmat(ROOT + 'data/rr_indices/baseline/' + filename)
                elif step == 'end ablation':
                    mat = loadmat(ROOT + 'data/rr_indices/end ablation/' + filename)
                    
                data_dict = mat['rr_indices_struct']
                feat_idx = list(data_dict).index(feat_name)
                data_matrix = np.stack(list(data_dict.values()))
                data_matrix = data_matrix.transpose()
                
                f = []
                p = []
                feat_names = data_dict.keys()
                p_norm = np.zeros((3, len(feat_names)))
                ss = ShuffleSplit(n_splits=3, test_size=0.2, random_state=0)
                classifier = LogisticRegression(penalty='l2', fit_intercept=True, random_state=0)
                
                tprs = []
                aucs = []
                mean_fpr = np.linspace(0, 1, 100)
                fig, ax = plt.subplots()
                
                for i, (train, test) in enumerate(ss.split(data_matrix)):
                    X_var = data_matrix[:,feat_idx].reshape(-1, 1)
                    classifier.fit(X_var[train], y[train])
                    viz = RocCurveDisplay.from_estimator(
                        classifier,
                        X_var[test],
                        y[test],
                        name="ROC fold {}".format(i),
                        pos_label=0,
                        alpha=0.3,
                        lw=1,
                        ax=ax,
                    )
                    plt.close()
                    
                    y_pred = classifier.predict_proba(X_var[test])
                    y_pred = y_pred[:,0]
                    fpr, tpr, thresholds = roc_curve(y[test], y_pred)
                    
                    interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
                    interp_thr = np.interp(mean_fpr, viz.fpr, thresholds)
                    interp_tpr[0] = 0.0
                    
                    J = interp_tpr - mean_fpr
                    idx = np.argmax(J)
                    opt_th = interp_thr[idx]
                    
                    tprs.append(interp_tpr)
                    aucs.append(viz.roc_auc)
                    
                mean_tpr = np.mean(tprs, axis=0)
                mean_tpr[-1] = 1.0
                mean_auc = auc(mean_fpr, mean_tpr)  #!!!!
                std_auc = np.std(aucs)
                
                J = mean_tpr - mean_fpr
                idx = np.argmax(J)
                opt_th = interp_thr[idx]
                sens = mean_tpr[idx]
                spec = 1 - mean_fpr[idx]

                df.append(
                {
                    'tau': tau,
                    'm': m,
                    'coef': coef,
                    'auc': mean_auc,
                    'std_auc': std_auc,
                    'sens': sens,
                    'spec': spec
                })
    
    df = pd.DataFrame(df)
    
    return df
                
    
    
def retrieve_boxplot_df(filt_type, norm_flag, var_to_plot, var_col_name, new_outcomes, outcome_names, pos_idx):
    sset_algos = ['seq','ecg12','ecg15']
    df = []

    # collect reconstruction errors for all subset selection algorithms, as well as for PCA with same number of PCs
    for n_sset_els in [8, 11]:
        for k, sset_algo in enumerate(sset_algos):
            # ecg12 and ecg15 consist of 8 and 11 electrodes, respectively (points on graph)
            if sset_algo == 'ecg12' and n_sset_els != 8:
                continue
            if sset_algo == 'ecg15' and n_sset_els != 11:
                continue

            # load in data
            ROOT = "C:/Users/amcca/switchdrive/PhD/frontiers_paper/"
            DATA_TYPE = 'complexity_indices_' + filt_type + '_'
            filename = DATA_TYPE + sset_algo + '_' + str(n_sset_els)
            if norm_flag:
                filename = filename  + '_normalized' 
            mat = loadmat(ROOT + 'data/complexity_indices/' + filename)
            data_dict = mat['complexity_indices_struct']

            # append subset algo and pca reconstruction error
            for pt_nb in pos_idx:  
                if sset_algo == 'seq':
                    if n_sset_els == 8:
                        sset_algo_name = '$SEQ_8$'
                    elif n_sset_els == 11:
                        sset_algo_name = '$SEQ_{11}$'
                elif sset_algo == 'ecg12':
                    sset_algo_name = '$ECG_8$'
                elif sset_algo == 'ecg15':
                    sset_algo_name = '$ECG_{11}$';

                df.append(
                {
                    'tau, m, coef': n_sset_els,
                    'Subset used': sset_algo_name,
                    var_col_name: data_dict[var_to_plot][pt_nb],
                    'Outcome': outcome_names[int(new_outcomes[pt_nb])]
                })


    # convert dictionary to pandas dataframe for use with seaborn
    df = pd.DataFrame(df)
    
    return df

def retrieve_boxplot_df_ndi(filt_type, norm_flag, var_to_plot, var_col_name, new_outcomes, outcome_names, pos_idx):
    df = []
    
    sset_algo = 'seq'
    n_sset_els = 8

    # load in data
    ROOT = "C:/Users/amcca/switchdrive/PhD/frontiers_paper/"
    DATA_TYPE = 'complexity_indices_' + filt_type + '_'
    filename = DATA_TYPE + sset_algo + '_' + str(n_sset_els)
    if norm_flag:
        filename = filename  + '_normalized' 
    mat = loadmat(ROOT + 'data/complexity_indices/' + filename)
    data_dict = mat['complexity_indices_struct']

    # append subset algo and pca reconstruction error
    for pt_nb in pos_idx:  
        df.append(
        {
            'x': 'x',
            var_col_name: data_dict[var_to_plot][pt_nb],
            'Outcome': outcome_names[int(new_outcomes[pt_nb])]
        })

    # convert dictionary to pandas dataframe for use with seaborn
    df = pd.DataFrame(df)
    
    return df


def retrieve_errcurve_df(filt_type, norm_flag, var_to_plot, pca_var, var_col_name, new_outcomes, outcome_names, pos_idx):
    sset_algos = ['seq','ecg12','ecg15','iter','random']
    sset_algo_names = ['$SEQ$', '$ECG_8$', '$ECG_{11}$', 'ITER', 'RAND']
    df = []
    
    # collect reconstruction errors for all subset selection algorithms, as well as for PCA with same number of PCs
    for n_sset_els in range(8,31):
        for k, sset_algo in enumerate(['seq','ecg12','ecg15','iter','random']):
            # ecg12 and ecg15 consist of 8 and 11 electrodes, respectively (points on graph)
            if sset_algo == 'ecg12' and n_sset_els != 8:
                continue
            if sset_algo == 'ecg15' and n_sset_els != 11:
                continue

            # load in data
            ROOT = "C:/Users/amcca/switchdrive/PhD/frontiers_paper/"
            DATA_TYPE = 'complexity_indices_' + filt_type + '_'
            filename = DATA_TYPE + sset_algo + '_' + str(n_sset_els)
            if norm_flag:
                filename = filename  + '_normalized' 
            mat = loadmat(ROOT + 'data/complexity_indices/' + filename)
            data_dict = mat['complexity_indices_struct']

            # append subset algo and pca reconstruction error
            for pt_nb in pos_idx:  
                df.append(
                {
                    'Number of PCs or electrodes in subset': n_sset_els,
                    'Subset algorithm': sset_algo_names[k],
                    var_col_name: data_dict[var_to_plot][pt_nb],
                    'Outcome': outcome_names[int(new_outcomes[pt_nb])]
                })


                if sset_algo == 'seq':
                    df.append(
                    {
                        'Number of PCs or electrodes in subset': n_sset_els,
                        'Subset algorithm': "$PCA$",
                        var_col_name: data_dict[pca_var][pt_nb],
                        'Outcome': outcome_names[int(new_outcomes[pt_nb])]
                    })


    # convert dictionary to pandas dataframe for use with seaborn
    df = pd.DataFrame(df)
    
    return df
    