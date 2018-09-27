#!/usr/bin/env python

import argparse
import numpy as np
import ROOT as r
import root_numpy as rnp
import sys
import yaml
import xgboost as xgb

from rootpy.plotting import Graph, Canvas
from sklearn import cross_validation
from sklearn.svm import LinearSVC
from sklearn import metrics 
from sklearn.externals import joblib

def parse_config(config_file) :

    print "Loading configuration from " + str(config_file)
    config = open(config_file, 'r')
    return yaml.load(config)

def main(): 

    # Parse all command line arguments using the argparse module
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-c", action='store', dest='config',
                        help="Configuration file.")
    args = parser.parse_args()

    # If a configuration file hasn't been specified, warn the user and exit.
    if not args.config :
        parser.error('A configuration file needs to be specified.')

    # Parse the configuration file
    config = parse_config(args.config)

    # Get the list of features that will be used to train the algorithm 
    if 'Features' not in config: 
        print '[ Trainer ] : A list of feature is required!'
        sys.exit(0)

    # Get the file containing the 'signal' and open it.  If a file isn't 
    # specified, warn the user and exit.
    if 'Signal' not in config: 
        print '[ Trainer ] : Signal file is required.'
        sys.exit(0)
    sig = rnp.root2array(config['Signal'], 'ecal_ntuple', branches=config['Features'])
    sig_trig = rnp.root2array(config['Signal'], 'ecal_ntuple', ['trigger_energy_sum'])
    sig_ntuple = rnp.root2array(config['Signal'], 'signal_ntuple')

    # Get the file containing the 'background' and open it.  If a file isn't 
    # specified, warn the user and exit.
    if 'Background' not in config: 
        print '[ Trainer ] : Background file is required.'
        sys.exit(0)

    bkg = rnp.root2array(config['Background'], 'ecal_ntuple', branches=config['Features'])
    bkg_trig = rnp.root2array(config['Background'], 'ecal_ntuple', ['trigger_energy_sum'])
    bkg_ntuple = rnp.root2array(config['Background'], 'pn_ntuple')

    if 'TriggerCut' not in config: 
        print '[ Trainer ]: An energy sum value is required by the trigger.'
        sys.exit(0)
    sig = sig[sig_trig['trigger_energy_sum'] < float(config['TriggerCut'][0])]
    bkg = bkg[bkg_trig['trigger_energy_sum'] < float(config['TriggerCut'][0])]
    sig_ntuple = sig_ntuple[sig_trig['trigger_energy_sum'] < float(config['TriggerCut'][0])]
    bkg_ntuple = bkg_ntuple[bkg_trig['trigger_energy_sum'] < float(config['TriggerCut'][0])]

    # Determine the number of events to train on.
    print '[ Trainer ]: Total signal: %s, Total background: %s' % (len(sig), len(bkg)) 
    t_events = min(len(sig), len(bkg))
    #t_events = 1000
    print '[ Trainer ] : Using %s events to train.' % t_events

    sig = rnp.rec2array(sig[:t_events])
    bkg = rnp.rec2array(bkg[:t_events])
    
    sig_en = np.arange(t_events)
    bkg_en = np.arange(t_events, 2*t_events)

    y_sig = np.ones(t_events)
    y_bkg = np.zeros(t_events)

    X = np.concatenate((sig, bkg))
    y = np.concatenate((y_sig, y_bkg))
    en = np.concatenate((sig_en, bkg_en))

    print en

    X_train, X_test, y_train, y_test, en_train, en_test = cross_validation.train_test_split(
            X, y, en, test_size=0.25, random_state=0)

    dtrain = xgb.DMatrix(X_train, label=y_train)
    dtest = xgb.DMatrix(X_test, label=y_test)

    param = {
        'objective': 'binary:logistic',
        'eta': 0.023,  # the training step for each iteration
        'max_depth': 10,  # the maximum depth of each tree
        'min_child_weight': 20,
        'silent': 1,  # logging mode - quiet
        'subsample':.9,
        'colsample_bytree':.85,
        'eval_metric':['auc', 'logloss', 'error'],
        'seed':1, 
        'nthread':3, 
        'verbosity':1, 
        'early_stopping_rounds':10
    }  
    
    #num_round = t_events  # the number of training iterations
    num_round = 1000  # the number of training iterations

    watchlist = [(dtest, 'eval')]
    evals_result = {}
    bst = xgb.train(param, dtrain, num_round, watchlist, verbose_eval=10, evals_result=evals_result)

    csv_arr = np.column_stack((evals_result['eval']['logloss'], 
                               evals_result['eval']['auc'],
                               evals_result['eval']['error']))
    np.savetxt(config['EvalResults_FileName'][0], csv_arr, delimiter=',')

    preds = bst.predict(dtest)

    preds_dic = {}
    for i in xrange(0, len(preds)): 
        preds_dic[en_test[i]] = preds[i]

    fpr, tpr, threshold = metrics.roc_curve(y_test, preds)

    roc_auc = metrics.auc(fpr, tpr)
    print 'Final Validation AUC = %s' % (roc_auc)

    csv_arr = np.column_stack((fpr, tpr))
    np.savetxt(config['ROC_FileName'][0], csv_arr, delimiter=',')

    joblib.dump(bst, config['PKL_FileName'][0], compress=True)

    sig_ntuple = sig_ntuple[:t_events]
    bkg_ntuple = bkg_ntuple[:t_events]

    sig_sel = []
    bkg_sel = []
    sig_pred = []
    bkg_pred = []
    for i in xrange(0, t_events): 
        if i in en_test: 
            sig_sel.append(True)
            sig_pred.append(preds_dic[i])
        else: sig_sel.append(False)

    sig_pred = np.array(sig_pred, dtype=[('pred', np.float64)])

    for i in xrange(t_events, 2*t_events): 
        if i in en_test: 
            bkg_sel.append(True)
            bkg_pred.append(preds_dic[i])
        else: bkg_sel.append(False)

    bkg_pred = np.array(bkg_pred, dtype=[('pred', np.float64)])
    
    sig_ntuple = sig_ntuple[np.array(sig_sel)]
    bkg_ntuple = bkg_ntuple[np.array(bkg_sel)]

    rnp.array2root(sig_ntuple, 'sig_test_ntuple.root', treename='signal_ntuple', mode='recreate')
    rnp.array2root(sig_pred, 'sig_test_preds.root', treename='signal_pred', mode='recreate')
    rnp.array2root(bkg_ntuple, 'bkg_test_ntuple.root', treename='pn_ntuple', mode='recreate')
    rnp.array2root(bkg_pred, 'bkg_test_preds.root', treename='pn_pred', mode='recreate')


if __name__ == "__main__":
    main()
