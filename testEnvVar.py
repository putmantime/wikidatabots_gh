__author__ = 'andra'
import os
import pprint

if '' in os.environ:
   print(os.environ['BLABLA'])

pprint.pprint(os.environ)
