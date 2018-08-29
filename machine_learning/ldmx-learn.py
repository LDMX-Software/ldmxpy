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
    sig_trig = rnp.root2array(config['Signal'], 'ecal_ntuple')
    #sig_trig = rnp.root2array(config['Signal'], 'trigger_ntuple')

    # Get the file containing the 'background' and open it.  If a file isn't 
    # specified, warn the user and exit.
    if 'Background' not in config: 
        print '[ Trainer ] : Background file is required.'
        sys.exit(0)
    bkg = rnp.root2array(config['Background'], 'ecal_ntuple', branches=config['Features'])
    bkg_trig = rnp.root2array(config['Background'], 'ecal_ntuple')
    #bkg_trig = rnp.root2array(config['Background'], 'trigger_ntuple')

    #sig = sig[sig_trig['triggered'] == 1]
    #bkg = bkg[bkg_trig['triggered'] == 1]

    sig = sig[sig_trig['trigger_energy_sum'] < 5650]
    bkg = bkg[bkg_trig['trigger_energy_sum'] < 5650]

    # Determine the number of events to train on.
    print '[ Trainer ]: Total signal: %s, Total background: %s' % (len(sig), len(bkg)) 
    t_events = min(len(sig), len(bkg))
    print '[ Trainer ] : Using %s events to train.' % t_events

    sig = rnp.rec2array(sig[:t_events])
    bkg = rnp.rec2array(bkg[:t_events])
    
    y_sig = np.ones(t_events)
    y_bkg = np.zeros(t_events)

    X = np.concatenate((sig, bkg))
    y = np.concatenate((y_sig, y_bkg))

    X_train, X_test, y_train, y_test = cross_validation.train_test_split(X, y, test_size=0.2, random_state=0)

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
        'eval_metric':'auc',
        'seed':1, 
        'nthread':1, 
        'verbosity':1, 
        'early_stopping_rounds':10
    }  
    
    num_round = 1000  # the number of training iterations

    bst = xgb.train(param, dtrain, num_round)

    preds = bst.predict(dtest)
    fpr, tpr, threshold = metrics.roc_curve(y_test, preds)

    roc_auc = metrics.auc(fpr, tpr)
    print 'Final Validation AUC = %s' % (roc_auc)

    csv_arr = np.column_stack((fpr, tpr))
    np.savetxt('roc_curve_2e_ecal_pn_ap_1mev.csv', csv_arr, delimiter=',')

    joblib.dump(bst, 'bst_model.pkl', compress=True)

if __name__ == "__main__":
    main()
