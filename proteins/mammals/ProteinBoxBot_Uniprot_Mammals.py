#!usr/bin/env python
# -*- coding: utf-8 -*-

'''
Author:Andra Waagmeester (andra@waagmeester.net)

This file is part of ProteinBoxBot.

ProteinBoxBot is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ProteinBoxBot is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ProteinBoxBot.  If not, see <http://www.gnu.org/licenses/>.
'''
# Load the path to the PBB_Core library
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../")
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../ProteinBoxBot_Core")

print (sys.path.append)

import ProteinBoxBot_Core.PBB_Core as PBB_Core

# Resource specific 
import protein
import traceback
from datetime import date, datetime, timedelta


try:
    if len(sys.argv) == 1:
        print("Please provide one of the following species as argument: human, mouse")
        print("Example: python ProteinBoxBot_Uniprot_Mammals.py human")
        sys.exit()

    proteome = protein.HumanProteome(sys.argv[1])


except Exception as err:
    print(traceback.format_exc())
    # client.captureException()  
