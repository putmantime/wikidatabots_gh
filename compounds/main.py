#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'Sebastian Burgstaller'
__licence__ = 'GPLv3'

"""
This file controls the drug bot. It can be run with two command line options:
User can choose to run either with aggregation of data or without. And the user can choose if the data file should be
generated from scratch or aggregation of data should just be continued.
"""

import sys
import argparse
from drug_data_aggregator import DrugDataAggregator
from drug_bot import DrugBot


def main():
    parser = argparse.ArgumentParser(description='Drug bot usage')
    parser.add_argument('--from-scratch', action='store_false', help='Build drug data file from scratch')
    parser.add_argument('--no-aggregation', action='store_true', help='Do not aggregate data, only run the bot')
    parser.add_argument('--user', action='store', help='Username on Wikidata', required=True)
    parser.add_argument('--pwd', action='store', help='Password on Wikidata', required=True)
    args = parser.parse_args()

    if args.no_aggregation:
        DrugBot(user=args.user, pwd=args.pwd)
        pass
    else:
        DrugDataAggregator(aggregate=args.from_scratch)
        DrugBot(user=args.user, pwd=args.pwd)

if __name__ == '__main__':
    sys.exit(main())
