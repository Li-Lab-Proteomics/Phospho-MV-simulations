import numpy as np
from sklearn.metrics import roc_curve, precision_recall_curve, auc, roc_auc_score
import pandas as pd

simple_implement=True
Nround=5 if simple_implement else 500
CD='../Phospho-MV-simulations-main'
sample_size=[5,10] if simple_implement==True else [5,10,30,50,80,100]
fclist=[1.5, 2] if simple_implement==True else [1.2, 1.5, 1.8, 2, 3, 5, 10]
zero_ratio=[0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
mnar_ratio=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
FDR=0.05


def roc_prc_scores(is_changed, p_vals):
    """
    Calculates AUROC, AUPRC, and standardized pAUROC scores
    :param is_changed: actual labels (vector of ones and zeros, n x 1)
    :param p_vals:     predicted labels (list of vectors of non-negative pvalues, smaller more significant, k x n)
    :return:           roc_auc, prc_aucs, pauc (all lists of floats)
    """
    roc_auc = []
    prc_auc = []
    pauc = []

    #for pval in p_vals.values.transpose():
    for pval in p_vals.transpose():
        # avoid trouble with -inf in roc_curve
        # replace 0 with an arbitrary small number
        pval = [p if p != 0 else 1e-100 for p in pval]

        predicted = - np.log(pval)

        # e.g. when values are missing because a test could not be applied
        if np.all(np.isnan(predicted)):
            # No valid p-vals
            pauc.append(np.nan)
            roc_auc.append(np.nan)
            prc_auc.append(np.nan)
            continue

        # Rank measurements with nan as lowest priority
        predicted[np.isnan(predicted)] = 0
        
        fpr, tpr, _ = roc_curve(is_changed, predicted)

        pauc_std=roc_auc_score(is_changed, predicted, max_fpr=FDR)

        pauc.append(pauc_std)
        roc_auc.append(auc(fpr, tpr))
        prec, rec, _ = precision_recall_curve(is_changed, predicted)
        prc_auc.append(auc(rec, prec))

    return roc_auc, prc_auc, pauc


def power_analysis(is_changed, pvals, alpha=0.05):
    """
    :param is_changed: actual labels (vector of ones and zeros, n x 1)
    :param p_vals:     predicted labels (list of vectors of non-negative pvalues, smaller more significant, n x k)
    :param alpha:      FWER, family-wise error rate, 0.05 by default in multipletests
    """
    fps = []
    tps = []
    Precision = []
    Recall = []
    for pval in pvals.transpose():
        #  e.g. when values are missing because the a test could not be applied
        if np.all(np.isnan(pval)):
            fp, tp, prec, rec = np.nan, np.nan, np.nan, np.nan

        else:
            # take NaN as negatives: (np.nan <= alpha) return False
            # False positives
            fp = np.sum(np.logical_not(is_changed) & (pval <= alpha))
            tp = np.sum((is_changed == 1) & (pval <= alpha))
            p = np.sum(is_changed)
            if p==len(is_changed):
                prec = np.nan
            elif fp+tp==0:
                prec = np.nan
            else:
                prec = tp/(tp+fp)
                
            rec = tp/p

        fps.append(fp)
        tps.append(tp)
        Precision.append(prec)
        Recall.append(rec)

    return fps, tps, Precision, Recall


Nmodel=17
method_name=['T-test','T-SampMin','T-bPCA','Wilcoxon','Wilcoxon-SampMin','Wilcoxon-bPCA','ModT','ModT-SampMin','ModT-bPCA','twoT','twoWilcox','SDA_robust','SDA','ZIG_2p','ZILN_2p','ZIG_DM','ZILN_DM']
result_all=pd.DataFrame()
for nsample in sample_size:
    for fc in fclist:
        for zr in zero_ratio:
            for mnar in mnar_ratio:
                for r in range(Nround):
                    r=r+1
                    # 17 different models
                    result_txt=CD+'/simu_results/round'+str(r)+'/Simu_'+str(nsample)+'_'+str(fc)+'_'+str(zr)+'_'+str(mnar)+'_for_'+str(r)+'_times.txt'
                    data=np.genfromtxt(result_txt,delimiter='\t',skip_header=1)
                    actual=data[:,data.shape[1]-1]
                    p_vals=data[:,0:(data.shape[1]-1)]
                    if np.all(actual==1):
                        auroc, auprc, pauroc=[np.nan]*Nmodel,[np.nan]*Nmodel,[np.nan]*Nmodel
                    else:
                        auroc, auprc, pauroc=roc_prc_scores(actual,p_vals)
                    FP, TP, precision, recall=power_analysis(actual, pvals=p_vals)
                    result_auc=np.array([auroc,pauroc,auprc,precision,recall])
                    result_scenario=np.array([nsample,fc,zr,mnar,r])
                    result_comb=np.hstack((np.tile(result_scenario,(Nmodel,1)),result_auc.T))                   
                    result_temp=pd.DataFrame(result_comb)
                    result_temp.insert(5,'model',method_name)
                    result_all=result_all.append(result_temp)
                    

result_all.columns=['Nsample','fc','zr','MNAR','round','model','AUROC','pAUROC','AUPRC','Precision','Recall']
output_dict={
    'Nsample':'int8',
    'fc':'float16',
    'zr':'float16',
    'MNAR':'float16',
    'round':'int16',
    'model':'object',
    'AUROC':'float64',
    'pAUROC':'float64',
    'AUPRC':'float64',
    'Precision':'float64',
    'Recall':'float64'}
result_all=result_all.astype(output_dict)
path_result=CD+'/Simu_results.txt'
result_all.to_csv(path_result,sep='\t',index=False)



