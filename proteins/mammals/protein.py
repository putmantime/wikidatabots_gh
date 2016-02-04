#!usr/bin/env python
# -*- coding: utf-8 -*-

"""
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
"""

__author__ = 'Andra Waagmeester and Jasper Koehorst'
__license__ = 'GPL'

import sys, os


sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../ProteinBoxBot_Core/")
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../")

print (sys.path)

import ProteinBoxBot_Core.PBB_Core as PBB_Core
import ProteinBoxBot_Core.PBB_Debug as PBB_Debug
import ProteinBoxBot_Core.PBB_login as PBB_login
import ProteinBoxBot_Core.PBB_settings as PBB_settings
from SPARQLWrapper import SPARQLWrapper, JSON
import requests
import traceback
from time import gmtime, strftime
import time
import pprint
import urllib.parse

global uniprot
uniprotURL = "http://sparql.uniprot.org/sparql/?query="

try:
    import simplejson as json
except ImportError as e:
    import json


def SPARQL(url, query, rformat):
    print ("Initialising SPARQL")
    prefixes = open("./resources/queries/prefixes").read()
    query = prefixes + query
    if rformat == "JSON":
        sparql = SPARQLWrapper(url)
        sparql.setQuery(prefixes+"\n"+query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        return results
    else:
        query = url + urllib.parse.quote_plus(query)
        result = requests.get(query + "&format="+rformat)
        results = result.json()
        return results

def GetWikidataIds(propertyId, wd_species):
        query = open("resources/queries/getIdInWikidata.sparql").read()
        query = query.replace("$1", propertyId)
        query = query.replace("$2", wd_species)
        return SPARQL("https://query.wikidata.org/bigdata/namespace/wdq/sparql", query, "JSON")

class DomainWikidataIDs(object):
    def __init__(self, object):
        """
        :rtype : basestring
        """
        if object == "human":
            self.taxon_id = 9606
            self.wd_taxon_id = "Q5"
            self.en_name = "human"
            self.fr_name = "humaine"
            self.nl_name = "menselijk"
            self.de_name = "Humanes"
        elif object == "mouse":
            self.taxon_id = 10090
            self.wd_taxon_id = "Q83310"
            self.en_name = "mouse"
            self.fr_name = "de souris"
            self.nl_name = "muizen"
            self.de_name = "Mause"
        elif object == "rat":
            self.taxon_id = 10116
            self.wd_taxon_id = "Q184224"
            self.en_name = "rat"
            self.fr_name = "de rat"
            self.nl_name = "rat"
            self.de_name = "Ratte"

class Proteins(object):
    def __init__(self, object):
        print("Getting all human proteins from Uniprot...")
        #Obtaining the right query
        query = open("resources/queries/proteinsByTaxon.sparql").read()
        query = query.replace("$1",object)
        #Perform the sparql query, if JSON sparqlwrapper is used otherwise, requests.get()
        prot_results = SPARQL(uniprotURL, query, "srj")

        self.uniprot_ids = []
        for protein in prot_results["results"]["bindings"]:
            item = dict()
            item["id"] = protein["protein"]["value"].replace("http://purl.uniprot.org/uniprot/", "")
            item["label"] = protein["protein_label"]["value"]
            self.uniprot_ids.append(item)

class Wikidata_Proteins(object):
    def __init__(self, object):
        print("Getting all proteins from Wikidata...")
        query = open("resources/queries/WikidataProteins")

class Uniprot_Record(object):
    def __init__(self, object):
        #Contains a $1 variable in the query which can be replaced by a uri / string / int etc of interest
        print (object["id"])
        query = open("resources/queries/uniprotAnnotations.sparql").read()
        query = query.replace("$1",object["id"])
        print (query)
        print ("Getting Uniprot records for a given protein")
        self.protein = SPARQL(uniprotURL, query,"srj")
        print(self.protein)
        if len(self.protein["results"]["bindings"]) == 0:
            raise Exception("Communication error on " + object)

class Go_Annotations_Uniprot(object):
    #Jasper: Not sure here... was (self):
    def __init__(self, object):
        query = open("resources/queries/GOTerms").read()
        query = query.replace("$1",object["id"])
        self.go_terms = SPARQL(uniprotURL,query,"srj")


class HumanProteome(object):
    def __init__(self, object):
        self.start = time.time()
        self.logincreds = PBB_login.WDLogin(PBB_settings.getWikiDataUser(), PBB_settings.getWikiDataPassword())
        self.domain_wd_ids = DomainWikidataIDs(object)
        
        print('Getting all "+object+" proteins from uniprot...')
        self.proteins = Proteins(object)

        print('Getting all "+object+" genes with a ncbi gene ID in Wikidata...')
        self.entrez_wikidata_ids = GetWikidataIds("P351", DomainWikidataIDs(object).wd_taxon_id)

        print('Getting all "+object+" proteins with a uniprot ID in Wikidata...')
        self.uniprot_wikidata_ids = GetWikidataIds("P352", DomainWikidataIDs(object).wd_taxon_id)


        InWikiData = PBB_Core.WDItemList(wdqQuery, wdprop="351")

        '''
        Below a mapping is created between entrez gene ids and wikidata identifiers.
        '''
        for geneItem in InWikiData.wditems["props"]["351"]:
            entrez_wikidata_ids[str(geneItem[2])] = geneItem[0]
        #Jasper: Proteins needs an object, is that self?

        for up in Proteins(self).uniprot_ids:
            print ("UP",up)
            try:
                protein_instance = Uniprot_Record(up)
                protein_instance["goTerms"] = Go_Annotations_Uniprot(up)
                protein_instance["logincreds"] = self.logincreds
                protein_instance["id"] = up["id"]
                protein_instance["start"] = self.start
                protein_instance["geneSymbols"] = genesymbolwdmapping
                protein_instance["entrez_wikidata_ids"] = entrez_wikidata_ids


            except Exception as e:
                print(traceback.format_exc())
                PBB_Core.WDItemEngine.log('ERROR',
                                          '{main_data_id}, "{exception_type}", "{message}", {wd_id}, {duration}'.format(
                                              main_data_id=up["id"],
                                              exception_type=type(e),
                                              message=e.__str__(),
                                              wd_id='-',
                                              duration=time.time() - self.start
                                          ))


class HumanProtein(object):
    def __init__(self, object):
        self.geneSymbols = object["geneSymbols"]
        self.logincreds = object["logincreds"]
        self.goTerms = object["goTerms"]
        self.version = object["results"]["bindings"][0]["upversion"]["value"]
        self.uniprot = object["results"]["bindings"][0]["uniprot"]["value"]
        self.uniprotId = object["id"]
        self.name = object["results"]["bindings"][0]["plabel"]["value"]
        self.start = object["start"]
        self.entrezWikidataIds = object["entrezWikidataIds"]

        if "gene_id" in object["results"]["bindings"][0].keys():
            self.gene_id = []
            for geneId in object["results"]["bindings"][0]["gene_id"]["value"].split(";"):
                if geneId != "":
                    self.gene_id.append(geneId)

        if "ecName" in object["results"]["bindings"][0].keys():
            self.ecname = []
            self.ecname.append(object["results"]["bindings"][0]["ecName"]["value"])

        self.alias = []
        for syn in object["results"]["bindings"][0]["upalias"]["value"].split(";"):
            if syn != "":
                self.alias.append(syn)
        if "pdbid" in object["results"]["bindings"][0].keys() and object["results"]["bindings"][0]["pdbid"][
            "value"] != "":
            self.pdb = []
            for pdbId in object["results"]["bindings"][0]["pdbid"]["value"].split(";"):
                self.pdb.append(pdbId.replace("http://rdf.wwpdb.org/pdb/", "").replace(" ", ""))
        if "refseqid" in object["results"]["bindings"][0].keys():
            self.refseq = []
            for refseqId in object["results"]["bindings"][0]["refseqid"]["value"].split(";"):
                self.refseq.append(refseqId.replace("http://purl.uniprot.org/refseq/", "").replace(" ", ""))
        if "ensemblp" in object["results"]["bindings"][0].keys() and object["results"]["bindings"][0]["ensemblp"][
            "value"] != "":
            self.ensemblp = []
            for ensP in object["results"]["bindings"][0]["ensemblp"]["value"].split(";"):
                self.ensemblp.append(ensP.replace("http://purl.uniprot.org/ensembl/", "").replace(" ", ""))

        # Prepare references
        refStatedIn = PBB_Core.WDItemID(value=2629752, prop_nr='P248', is_reference=True)
        refStatedIn.overwrite_references = True
        refURL = "http://www.uniprot.org/uniprot/" + self.uniprotId + ".txt?version=" + str(self.version)
        refReferenceURL = PBB_Core.WDUrl(value=refURL, prop_nr='P854', is_reference=True)
        refReferenceURL.overwrite_references = True
        refImported = PBB_Core.WDItemID(value=905695, prop_nr='P143', is_reference=True)
        refImported.overwrite_references = True
        timeStringNow = strftime("+%Y-%m-%dT00:00:00Z", gmtime())
        refRetrieved = PBB_Core.WDTime(timeStringNow, prop_nr='P813', is_reference=True)
        refRetrieved.overwrite_references = True
        protein_reference = [[refStatedIn, refImported, refRetrieved, refReferenceURL]]

        references = dict()
        proteinPrep = dict()

        # P279 = subclass of
        proteinPrep['P279'] = [PBB_Core.WDItemID(value="Q8054", prop_nr='P279', references=protein_reference)]

        # P703 = found in taxon
        proteinPrep['P703'] = [PBB_Core.WDItemID(value="Q5", prop_nr='P703', references=protein_reference)]

        # P352 = UniprotID
        proteinPrep['P352'] = [PBB_Core.WDString(value=self.uniprotId, prop_nr='P352', references=protein_reference)]

        # P591 = ec number
        if "ecname" in vars(self):
            proteinPrep['P591'] = []
            for i in range(len(self.ecname)):
                proteinPrep['P591'].append(
                    PBB_Core.WDString(value=self.ecname[i], prop_nr='P591', references=protein_reference))

        # P638 = PDBID
        if "pdb" in vars(self) and len(self.pdb) > 0:
            proteinPrep['P638'] = []
            for i in range(len(self.pdb)):
                proteinPrep['P638'].append(
                    PBB_Core.WDString(value=self.pdb[i], prop_nr='P638', references=protein_reference))

        # P637 = Refseq Protein ID
        if "refseq" in vars(self) and len(self.refseq) > 0:
            proteinPrep['P637'] = []
            for i in range(len(self.refseq)):
                proteinPrep['P637'].append(
                    PBB_Core.WDString(value=self.refseq[i], prop_nr='P637', references=protein_reference))

        # P705 = Ensembl Protein ID
        if "ensemblp" in vars(self) and len(self.ensemblp) > 0:
            proteinPrep['P705'] = []
            for i in range(len(self.ensemblp)):
                proteinPrep['P705'].append(
                    PBB_Core.WDString(value=self.ensemblp[i], prop_nr='P705', references=protein_reference))

        """
        # P686 = Gene Ontology ID
        proteinPrep["P680"] = []
        proteinPrep["P681"] = []
        proteinPrep["P682"] = []

        for result in self.goTerms["results"]["bindings"]:

            statement = [
                    PBB_Core.WDString(value=result["go"]["value"].replace("http://purl.obolibrary.org/obo/GO_", "GO:"),
                                      prop_nr='P686', references=protein_reference)]
            goWdPage = PBB_Core.WDItemEngine(item_name=result["goLabel"]["value"], data=statement,
                                                 server="www.wikidata.org", domain="proteins")
            if goWdPage.get_description() == "":
                goWdPage.set_description("Gene Ontology term")
            js = goWdPage.get_wd_json_representation()
            goWdId = goWdPage.write(self.logincreds)

            if result["parentLabel"]["value"] == "molecular_function":
                exists = False
                for i in range(len(proteinPrep["P680"])):
                    if proteinPrep["P680"][i].value == goWdId:
                        exists = True
                if not exists:
                    proteinPrep["P680"].append(
                        PBB_Core.WDItemID(value=goWdId, prop_nr='P680', references=protein_reference))
            if result["parentLabel"]["value"] == "cellular_component":
                exists = False
                for i in range(len(proteinPrep["P681"])):
                    if proteinPrep["P681"][i].value == goWdId:
                        exists = True
                if not exists:
                    proteinPrep["P681"].append(
                        PBB_Core.WDItemID(value=goWdId, prop_nr='P681', references=protein_reference))
            if result["parentLabel"]["value"] == "biological_process":
                exists = False
                for i in range(len(proteinPrep["P682"])):
                    if proteinPrep["P682"][i].value == goWdId:
                        exists = True
                if not exists:
                    proteinPrep["P682"].append(
                        PBB_Core.WDItemID(value=goWdId, prop_nr='P682', references=protein_reference))
        """

        # P702 = Encoded by
        if "gene_id" in vars(self) and len(self.gene_id) > 0:
            proteinPrep['P702'] = []
            proteinPrep['P702'].append(
                PBB_Core.WDItemID(value=self.entrezWikidataIds[
                    self.gene_id[0].replace("http://purl.uniprot.org/geneid/", "").replace(" ", "")], prop_nr='P702',
                                  references=protein_reference))

        proteinData2Add = []
        for key in proteinPrep.keys():
            for statement in proteinPrep[key]:
                proteinData2Add.append(statement)
                print(statement.prop_nr, statement.value)
        wdProteinpage = PBB_Core.WDItemEngine(item_name=self.name, data=proteinData2Add, server="www.wikidata.org",
                                                  domain="proteins", append_value=['P279'])

        if len(self.alias) > 0:
            wdProteinpage.set_aliases(aliases=self.alias, lang='en', append=True)
        if wdProteinpage.get_description() == "":
            wdProteinpage.set_description(description='human protein', lang='en')
        if wdProteinpage.get_description(lang="de") == "":
            wdProteinpage.set_description(description='humanes Protein', lang='de')
        if wdProteinpage.get_description(lang="nl") == "":
            wdProteinpage.set_description(description='menselijk eiwit', lang='nl')
        if wdProteinpage.get_description(lang="fr") == "" or wdProteinpage.get_description(lang="fr") == "protéine":
            wdProteinpage.set_description(description='protéine humaine', lang='fr')

        self.wd_json_representation = wdProteinpage.get_wd_json_representation()
        # PBB_Debug.prettyPrint(self.wd_json_representation)
        #wdProteinpage.write(self.logincreds)
        print(wdProteinpage.wd_item_id)

        PBB_Core.WDItemEngine.log('INFO', '{main_data_id}, "{exception_type}", "{message}", {wd_id}, {duration}'.format(
            main_data_id=self.uniprotId,
            exception_type='',
            message=f.name,
            wd_id=self.wdid,
            duration=time.time() - self.start
        ))
        print("===============")
