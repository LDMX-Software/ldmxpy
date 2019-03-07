#!/usr/bin/env python

import argparse
import importlib
import Event as e
import ROOT as r
import os
import sys
import yaml

from rootpy.io import root_open 

def parse_config(config_file) :

    print "Loading configuration from " + str(config_file)
    config = open(config_file, 'r')
    return yaml.load(config)

def main() : 
   
    # Parse all command line arguments using the argparse module
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-c", action='store', dest='config',
                        help="Configuration file.")
    parser.add_argument("-n", action='store', dest='n_events', 
                        help="Total number of events.")
    parser.add_argument('-p', action='store', dest='n_print', 
                         help='Freqency of event number printing.')
    args = parser.parse_args()

    if not args.config :
        parser.error('A configuration file needs to be specified.')

    n_events = 0
    if args.n_events: n_events = args.n_events

    n_print = 1000
    if args.n_print: n_print = int(args.n_print) 

    # Parse the configuration file
    config = parse_config(args.config)
    
    analyses = config["Analyses"]
    analyses_instances = []
    for analysis in analyses : 
        analysis_module_name, analysis_class_name = analysis.rsplit(".", 1)
        print "[ ldmxpy ]: Adding analysis ==> Module: %s Class: %s" % (analysis_module_name, analysis_class_name)
        analysis_class = getattr(importlib.import_module(analysis_module_name), analysis_class_name)
        analyses_instances.append(analysis_class())

    files = []
    if "Files" in config:
        files = config["Files"]
    elif "FileList" in config: 
        flist = open(config["FileList"][0], 'r')
        for f in flist: 
            files.append(f.strip())

    ofile_path = 'analysis.root'
    if 'OutputFile' in config: 
        ofile_path = config['OutputFile'][0]

    ofile = root_open(ofile_path, 'recreate')
    
    params = {}
    if 'Parameters' in config: 
        params = config['Parameters']

    for analyses in analyses_instances : 
        analyses.initialize(params)
   
    tree_name = 'LDMX_Events'
    if 'TreeName' in config: 
        tree_name = config['TreeName'][0]

    event = e.Event(config)
    # Loop through all of the ROOT files and process them.
    for rfile_path in files :
        print 'Processing file %s' % rfile_path
        event.load_file(rfile_path, tree_name)
        
        event_counter = 0
        while event.next_event():
            for analysis in analyses_instances:
                analysis.process(event)
            event_counter += 1
            
            if event_counter == int(n_events): 
                print 'Hit event limit'
                break

            if event_counter%n_print == 0: 
                print '[ ldmxpy ]: >> Event >> %s >>' % event_counter

        print "Total number of events processed: %s" % event_counter
        event.close_file()

    for analyses in analyses_instances : 
        analyses.finalize()
        
    ofile.close()

if __name__ == "__main__":
    main()
