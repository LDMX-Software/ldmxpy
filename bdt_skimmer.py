#!/usr/bin/python

import argparse
import importlib
import Event as e
import ROOT as r
import os
import sys

def parse_config(config_file) :

    print "Loading configuration from " + str(config_file)
    config = open(config_file, 'r')
    return yaml.load(config)

def main() : 
   
    # Parse all command line arguments using the argparse module
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input_file", help="LDMX recon file")
    parser.add_argument("-l", "--ldmx_lib",  help="Location of LDMX-SW library")
    args = parser.parse_args()

    if not args.input_file:
        parser.error('A list of files to skim needs to be specified.')
   
    # Get the path for the event lib
    if not args.ldmx_lib:
        parser.error('The location of the LDXM-SW event library was not specified.')
    ldmx_lib_path = args.ldmx_lib

    event = e.Event(ldmx_lib_path)

    print 'Processing file %s' % args.input_file
    output_path = args.input_file.strip()[:-5]
    output_path = output_path[output_path.rfind('/') + 1:]
    output_path += '_skim.root'
    print 'Output Path: %s' % output_path 
    event.load_file(args.input_file.strip())
    rfile_new = r.TFile(output_path, 'recreate')
    tree_new = event.get_tree().CloneTree(0)
    while event.next_event():
        ecal_hits = event.get_collection('EcalSimHits')
        energy_counter = 0
        for ecal_hit in ecal_hits: 
            energy_counter += ecal_hit.getEdep()

        hcal_hits = event.get_collection('HcalHits')
        pe_counter = 0
        for hcal_hit in hcal_hits: 
            pe_counter += hcal_hit.getPE()

        if (energy_counter < 25) & (pe_counter < 24): tree_new.Fill()

    tree_new.AutoSave()
    event.close_file()
    rfile_new.Close()

if __name__ == "__main__":
    main() 
