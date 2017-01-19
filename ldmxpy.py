#!/usr/bin/python

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
    args = parser.parse_args()

    if not args.config :
        parser.error('A configuration file needs to be specified.')

    
    #ldmx_lib_path = os.environ['LDMX_SW_DIR'] + "/install/lib/libEvent.so"

    # Parse the configuration file
    config = parse_config(args.config)
    
    # Get the path for the event lib
    ldmx_lib_path = config['EventLib'][0]

    analyses = config["Analyses"]
    analyses_instances = []
    for analysis in analyses : 
        analysis_module_name, analysis_class_name = analysis.rsplit(".", 1)
        print "[ ldmxpy ]: Adding analysis ==> Module: %s Class: %s" % (analysis_module_name, analysis_class_name)
        analysis_class = getattr(importlib.import_module(analysis_module_name), analysis_class_name)
        analyses_instances.append(analysis_class())

    event = e.Event(ldmx_lib_path)
    # Loop through all of the ROOT files and process them.
    for rfile_path in config["Files"] :
        print 'Processing file %s' % rfile_path
        event.load_file(rfile_path)

        while event.next_event():
            for analysis in analyses_instances:
                analysis.process(event)
            
        event.close_file()
    for analyses in analyses_instances : 
        analyses.finalize()

if __name__ == "__main__":
    main()
