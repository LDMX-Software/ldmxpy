#!/usr/bin/env python

import argparse
import logging
import uproot4
import yaml

import numpy as np
import pandas as pd

logging.basicConfig(format='[ ldmxpy ][ %(levelname)s ]: %(message)s',
                    level=logging.INFO)

def main():
 
  # Parse all command line arguments using the argparse module
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('-c', action='store', dest='config',
                      help='Configuration in yaml format.')
  parser.add_argument('-f', action='store', dest='file', 
                      help='The input file to process.')
  args = parser.parse_args()
    
  # If a configuration file was not specified, warn the user and exit.
  if not args.config :
    parser.error('A configuration file needs to be specified.')
  
  if not args.file: 
    parser.error('A ROOT file to process needs to be specified.')

  # Parse the configuration file.
  logging.info('Parsing configuration located at %s' % args.config.strip()) 
  config = yaml.load(open(args.config.strip(), 'r'), Loader=yaml.FullLoader)
 
  logging.info('Processing file %s' % args.file.strip())
  tree = uproot4.open('%s:LDMX_Events' % args.file.strip())

  data_dict = {}
  for key, values in config['variables'].items():
    for value in values:
      a = tree[key][value].array(library='np')
      if a[0].size != 1:
        for i in range(0, a[0].size):
          tmp = []
          for j in range(0, a.size):
            tmp.append(a[j][i])
          data_dict['%s_%s' % (value[:-1], i)] = np.array(tmp)
      else:
        data_dict[value[:-1]] = tree[key][value].array(library='np')

  print(config['output'])
  data_frame = pd.DataFrame.from_dict(data_dict)
  data_frame.to_csv(config['output'].strip(), sep=',', index=False)

if __name__ == "__main__":
    main()
