__author__ = 'andra'
import os
import pprint

if 'WD_API' in os.environ:
   print(os.environ['WD_API'])

pprint.pprint(os.environ)
