#!usr/bin/env python
# -*- coding: utf-8 -*-

'''
Authors: 
  Andra Waagmeester (andra' at ' micelio.be)

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

__author__ = 'Andra Waagmeester'
__license__ = 'GPL'

import sys
sys.path.append("/Users/andra/wikidatabots/ProteinBoxBot_Core")
import PBB_Core
import PBB_Debug
import PBB_login
import PBB_settings
import urllib
import urllib3
import certifi

import sys
import mygene_info_settings

try:
    import simplejson as json
except ImportError as e:
    import json
    
class human_genome():
    def __init__(self):
        self.content = json.loads(self.download_human_genes())
        self.gene_count = self.content["total"]
        self.genes = self.content["hits"]
        self.logincreds = PBB_login.WDLogin(PBB_settings.getWikiDataUser(), PBB_settings.getWikiDataPassword())
        
        entrezWikidataId = dict()
        print "Getting all entrez genes in Wikidata"
        InWikiData = PBB_Core.WDItemList("CLAIM[703:5] AND CLAIM[351]", "351")
        for geneItem in InWikiData.wditems["props"]["351"]:
            entrezWikidataId[str(geneItem[2])] = geneItem[0]
        for gene in self.genes:
            geneClass = human_gene(gene)
            if str(geneClass.entrezgene) in entrezWikidataId.keys():
                geneClass.wdid = entrezWikidataId[str(geneClass.entrezgene)]
            else: 
                geneClass.wdid = None  

    def download_human_genes(self):
        """
        Downloads the latest list of human genes from mygene.info through the URL specified in mygene_info_settings
        """
        # request = urllib2.Request(mygene_info_settings.getHumanGenesUrl())
        urllib.urlretrieve (mygene_info_settings.getHumanGenesUrl(), "human_genes.json")
        file = open("human_genes.json", 'r')
        return file.read()
        
class human_gene(object):
    def __init__(self, object):
        self.content = object
        self.entrezgene = object["entrezgene"]
        self.name = object["name"]
        self.symbol = object["symbol"]
        self.wdid = ''
        gene_annotations = json.loads(self.annotate_gene())
        PBB_Debug.prettyPrint(gene_annotations)
        print gene_annotations["_id"]
        self.annotationstimestamp = gene_annotations["_timestamp"]
        if "alias" in gene_annotations.keys(): 
            self.synonyms = gene_annotations["alias"]
        else:
            self.synonyms = None
        self.ensembl_transcript = None
        self.ensembl_gene = None
        if self.wdid != None:
            wdPage = PBB_Core.WDItemEngine(self.wdid, self.name)
        else:
            wdPage = PBB_Core.WDItemEngine('', self.name)
        print "Test:"
        PBB_Debug.prettyPrint(wdPage.wd_json_representation)
        print self.wdid
        sys.exit()
        ''' if "ensembl" in gene_annotations.keys():
            if "gene" in gene_annotations["ensembl"].keys():
                self.ensembl_gene = gene_annotations["ensembl"]["gene"]
            if "transcript" in gene_annotations["ensembl"].keys():
                self.ensembl_transcript = gene_annotations["ensembl"]["transcript"]
        '''        
        
    def annotate_gene(self):
        "Get gene annotations from mygene.info"
        http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED', ca_certs=certifi.where())
        request = http.request("GET", mygene_info_settings.getGeneAnnotationsURL()+str(self.entrezgene))
        return request.data
        

        
        
        
        
        