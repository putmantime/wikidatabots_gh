__author__ = 'timputman'
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../ProteinBoxBot_Core")
from SPARQLWrapper import SPARQLWrapper, JSON
import csv
import requests
import urllib.request
import pprint
import PBB_Core
import PBB_login
from time import gmtime, strftime
import copy

__author__ = 'timputman'



class WDProp2QID_SPARQL(object):
    def __init__(self, prop='', string=''):
        self.qid = self.SPARQL_for_qidbyprop(prop, string)

    def SPARQL_for_qidbyprop(self, prop, string):
        """
        :param prop: 'P351' Entrez gene id (ex. print( SPARQL_for_qidbyprop('P351','899959')))
        :param string: '899959' String value
        :return: QID Q21514037
        """
        sparql = SPARQLWrapper("https://query.wikidata.org/bigdata/namespace/wdq/sparql")
        prefix = 'PREFIX wdt: <http://www.wikidata.org/prop/direct/>'
        arguments = '?gene wdt:{} "{}"'.format(prop, string)
        select_where = 'SELECT * WHERE {{{}}}'.format(arguments)
        query = prefix + " " + select_where
        sparql.setQuery(query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        final_qid =[]
        try:
            rawqid = results['results']['bindings'][0]['gene']['value']
            qid_list = rawqid.split('/')
            final_qid.append(qid_list[-1])
        except Exception:
            final_qid.append('None')
        return final_qid[0]


class WDQID2Label_SPARQL(object):
    def __init__(self, qid=''):
        self.label = self.qid2label(qid)

    def qid2label(self, qid):
        sparql = SPARQLWrapper("https://query.wikidata.org/bigdata/namespace/wdq/sparql")
        prefix = 'PREFIX wd: <http://www.wikidata.org/entity/>'
        arguments = ' wd:{} rdfs:label ?label. Filter (LANG(?label) = "en") .'.format(qid)
        select_where = 'SELECT ?label WHERE {{{}}}'.format(arguments)
        query = prefix + " " + select_where
        sparql.setQuery(query)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        final_qid =[]
        try:
            rawqid = results['results']['bindings'][0]['label']['value']
            final_qid.append(rawqid)
        except Exception:
            final_qid.append('None')
        return final_qid[0]


class NCBIReferenceGenomes(object):
    def __init__(self):
        self.tid_list = self.get_ref_microbe_taxids()

    def get_ref_microbe_taxids(self):
        """
        Download the latest bacterial genome assembly summary from the NCBI genome ftp site
        and generate a list of relevant data for strain items based on taxids of the bacterial reference genomes.

        :return:
        """
        assembly = urllib.request.urlopen("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt")
        datareader = csv.reader(assembly.read().decode().splitlines(), delimiter="\t")

        ref_org_list = []
        for row in datareader:
            if row[4] == 'reference genome':
                org = row[7].split()
                strain = " ".join(org[2:])
                genus = org[0]
                species = org[1]
                species_tid = row[6]
                strain_tid = row[5]
                organism_data = [species_tid, strain_tid, genus, species, strain]
                ref_org_list.append(organism_data)
        return ref_org_list


class StrainDataParser(object):
    def __init__(self,strain_data):
        self.strain_data = strain_data
        self.species_taxid = str(self.strain_data[0])
        self.strain_taxid = self.strain_data[1]
        self.taxon_name = " ".join(self.strain_data[2:])
        self.genus = self.strain_data[2]
        self.species = self.strain_data[3]
        self.strain = self.strain_data[4]


class GeneDataParser(object):
    def __init__(self, gene_data):
        self.gene_data = gene_data
        self.RS_genome = self.gene_data['RSgenomic']
        self.strand = str(self.gene_data['strand'])
        self.type_of_gene = self.gene_data['type_of_gene']
        self.name = self.gene_data['name']
        self.gene_symbol = self.gene_data['gene_symbol']
        self.geneid = str(self.gene_data['_geneid'])
        self.genstop = str(self.gene_data['genstop'])
        self.genstart = str(self.gene_data['genstart'])
        self.protein_symbol = self.gene_data['protein_symbol']
        self.RS_protein = self.gene_data['RSprotein']
        self.uniprotid = self.gene_data['uniprot']
        self.GO_BP = self.gene_data['biological_process']
        self.GO_MF = self.gene_data['molecular_function']
        self.GO_CC = self.gene_data['cell_component']
        self.ec_number = str(self.gene_data['ec_number'])


class MyGeneInfoRestBatchQuery(object):
    def __init__(self, strain_data):
        self.strain_taxid = StrainDataParser(strain_data).strain_taxid
        self.gene_record = self.mygeneinfo_download()

    def mygeneinfo_download(self):
        url = 'http://mygene.info/v2/query/'
        params = dict(q="__all__", species=self.strain_taxid, entrezonly="true", size="20", fields="all")
        r = requests.get(url=url, params=params)
        hits = r.json()

        mgi_data = []

        for i in hits['hits']:
            if i['type_of_gene'] != 'protein-coding':
                continue
            else:
                uniprot = ''
                if 'uniprot' in i.keys():
                    if 'Swiss-Prot' in i['uniprot']:
                        if isinstance(i['uniprot']['Swiss-Prot'],str):
                            uniprot = i['uniprot']['Swiss-Prot']
                        else:
                            uniprot = i['uniprot']['Swiss-Prot'][0]

                    if 'TrEMBL' in i['uniprot']:
                        if isinstance(i['uniprot']['TrEMBL'],str):
                            uniprot = i['uniprot']['TrEMBL']
                        else:
                            uniprot = i['uniprot']['TrEMBL'][0]
                else:
                    continue

                wd_data = {
                    '_geneid': str(i['entrezgene']),
                    'type_of_gene': i['type_of_gene'],
                    'gene_symbol': i['symbol'],
                    'protein_symbol': i['symbol'].upper(),
                    'name': i['name'],
                    'uniprot': uniprot,
                    'RSprotein': i['refseq']['protein'],
                    'genstart': str(i['genomic_pos']['start']),
                    'genstop': str(i['genomic_pos']['end']),
                    'strand': i['genomic_pos']['strand'],
                    'RSgenomic': i['genomic_pos']['chr'],
                    'biological_process': '',
                    'cell_component': '',
                    'molecular_function': '',
                    'ec_number': ''
                }
                mgi_data.append(wd_data)
        return mgi_data


class UniProtRESTBatchQuery(object):
    def __init__(self, strain_data):
        self.strain_taxid = StrainDataParser(strain_data).strain_taxid
        self.enzyme_record = self.download_taxon_protein_GO()
    def download_taxon_protein_GO(self):
        """
        Downloads the latest list of microbial proteins from uniprot, taxon specified by the strain takid provided.
        """

        url = 'http://www.uniprot.org/uniprot/'

        params = dict(query=('organism:' + self.strain_taxid), format='tab',
                      columns='id,go(biological process),go(cellular component),go(molecular function),ec')
        r = requests.get(url=url, params=params)
        go_terms = r.text
        datareader = csv.reader(go_terms.splitlines(), delimiter="\t")
        uniprot_data = []

        for i in datareader:
            go_dict = {
                'uniprot': i[0],
                'biological_process': i[1], # biological process P682
                'cell_component': i[2], # cell component P681
                'molecular_function': i[3], # molecular function P680
                'ec_number': i[4]
            }
            uniprot_data.append(go_dict)
        return uniprot_data


class MGI_UNIP_Merger(object):
    def __init__(self, mgi, unip):
        self.mgi_data = mgi
        self.unip_data = unip
        self.mgi_unip_dict = self.combine_mgi_uniprot_dicts()

    def combine_mgi_uniprot_dicts(self):
        all_data = []
        for m in self.mgi_data:
            for u in self.unip_data:
                if m['uniprot'] == u['uniprot']:
                    m.update(u)
            all_data.append(m)
        return all_data


class WDGeneItem(object):
    def __init__(self, gene_record, strain_record):
        self.gene_record = GeneDataParser(gene_record)
        self.strain_record = StrainDataParser(strain_record)
        self.final_statements = []

    def gene_item(self):
        """
        Write gene items using dictionary from wd_parse_mgi(self):
        :return:
        """
        item_name = self.gene_record.name + self.gene_record.gene_symbol
        alias_list = self.gene_record.gene_symbol
        strain_qid = WDProp2QID_SPARQL(prop='P685', string=self.strain_record.strain_taxid).qid
        strain_label = WDQID2Label_SPARQL(qid=strain_qid).label
        description = "microbial gene found in " + strain_label

        NCBIgenerefStated = PBB_Core.WDItemID(value='Q20641742', prop_nr='P248', is_reference=True)
        NCBIgenerefStated.overwrite_references = True
        timeStringNow = strftime("+%Y-%m-%dT00:00:00Z", gmtime())
        refRetrieved = PBB_Core.WDTime(str(timeStringNow), prop_nr='P813', is_reference=True)
        refRetrieved.overwrite_references = True
        NCIB_gene_reference = [NCBIgenerefStated, refRetrieved]

        statements = dict()
        statements['found_in'] = [PBB_Core.WDItemID(value=strain_qid, prop_nr='P703', references=[copy.deepcopy(NCIB_gene_reference)])] #Found in taxon
        statements['geneid'] = [PBB_Core.WDString(value=self.gene_record.geneid, prop_nr='P351', references=[copy.deepcopy(NCIB_gene_reference)])]
        statements['subclass'] = [PBB_Core.WDItemID(value='Q7187', prop_nr='P279', references=[copy.deepcopy(NCIB_gene_reference)])]
        statements['genomic_start'] = [PBB_Core.WDString(value=self.gene_record.genstart, prop_nr='P644', references=[copy.deepcopy(NCIB_gene_reference)])]
        statements['genomic_stop'] = [PBB_Core.WDString(value=self.gene_record.genstop, prop_nr='P645', references=[copy.deepcopy(NCIB_gene_reference)])]

        for key in statements.keys():
            for statement in statements[key]:
                self.final_statements.append(statement)
        wd_item_gene = PBB_Core.WDItemEngine(item_name=item_name, domain='genes', data=self.final_statements)
        wd_item_gene.set_label(item_name)
        wd_item_gene.set_description(description)
        wd_item_gene.set_aliases(alias_list)
        pprint.pprint(wd_item_gene.get_wd_json_representation())
        #wd_item_gene.write(login)


class WDProteinItem(object):
    def __init__(self, gene_record, strain_record):
        self.gene_record = GeneDataParser(gene_record)
        self.strain_record = StrainDataParser(strain_record)
        self.final_statements = []

    def protein_item(self):

        item_name = self.gene_record.name + self.gene_record.protein_symbol
        alias_list = self.gene_record.protein_symbol
        strain_qid = WDProp2QID_SPARQL(prop='P685', string=self.strain_record.strain_taxid).qid
        strain_label = WDQID2Label_SPARQL(qid=strain_qid).label
        description = "microbial protein found in " + strain_label

        uniprotrefStated = PBB_Core.WDItemID(value='Q905695', prop_nr='P248', is_reference=True)
        uniprotrefStated.overwrite_references = True
        timeStringNow = strftime("+%Y-%m-%dT00:00:00Z", gmtime())
        refRetrieved = PBB_Core.WDTime(str(timeStringNow), prop_nr='P813', is_reference=True)
        refRetrieved.overwrite_references = True
        uniprot_protein_reference = [uniprotrefStated, refRetrieved]


        statements = dict()
        statements['found_in'] = [PBB_Core.WDItemID(value=strain_qid, prop_nr='P703', references=[copy.deepcopy(uniprot_protein_reference)])] #Found in taxon
        statements['uniprot'] = [PBB_Core.WDString(value=self.gene_record.uniprotid, prop_nr='P352', references=[copy.deepcopy(uniprot_protein_reference)])]
        statements['proteinid'] = [PBB_Core.WDString(value=self.gene_record.RS_protein, prop_nr='P637', references=[copy.deepcopy(uniprot_protein_reference)])]
        statements['subclass'] = [PBB_Core.WDItemID(value='Q8054', prop_nr='P279', references=[copy.deepcopy(uniprot_protein_reference)])]
        statements['molecular_function'] = []
        statements['cell_component'] = []
        statements['biological_process'] = []

        if self.gene_record.ec_number:
            statements['ec_number'] = [PBB_Core.WDString(value=self.gene_record.ec_number, prop_nr='P591', references=[copy.deepcopy(uniprot_protein_reference)])]

        if self.gene_record.GO_MF:
            mflist = self.gene_record.GO_MF.split(";")  # P680
            mfgoqidlist= []
            for item in mflist:
                mfdict = {}
                mf = item.split()
                mfdict['mf_goid'] = mf[-1][1:-1]
                mfdict['mf_goterm'] = " ".join(mf[:-1])
                mfdict['mf_qoqid'] =  WDProp2QID_SPARQL(prop='P686', string= mfdict['mf_goid']).qid
                mfdict['mf_golabel'] = WDQID2Label_SPARQL(qid=mfdict['mf_qoqid']).label

                if mfdict['mf_goterm'] == mfdict['mf_golabel']:
                    mfgoqidlist.append(mfdict['mf_qoqid'])
                    #add log statement
            statements['molecular_function'].append(PBB_Core.WDItemID(value=mfgoqidlist, prop_nr='P680', references=[copy.deepcopy(uniprot_protein_reference)]))

        if self.gene_record.GO_CC:
            cclist = self.gene_record.GO_CC.split(";")  # P680
            ccgoqidlist= []
            for item in cclist:
                ccdict = {}
                cc = item.split()
                ccdict['cc_goid'] = cc[-1][1:-1]
                ccdict['cc_goterm'] = " ".join(cc[:-1])
                ccdict['cc_qoqid'] =  WDProp2QID_SPARQL(prop='P686', string= ccdict['cc_goid']).qid
                ccdict['cc_golabel'] = WDQID2Label_SPARQL(qid=ccdict['cc_qoqid']).label

                if ccdict['cc_goterm'] == ccdict['cc_golabel']:
                    ccgoqidlist.append(ccdict['cc_qoqid'])
                      #add log statement
            statements['cell_component'].append(PBB_Core.WDItemID(value=ccgoqidlist, prop_nr='P680', references=[copy.deepcopy(uniprot_protein_reference)]))

        if self.gene_record.GO_BP:
            bplist = self.gene_record.GO_BP.split(";")  # P680
            bpgoqidlist= []
            for item in bplist:
                bpdict = {}
                bp = item.split()
                bpdict['bp_goid'] = bp[-1][1:-1]
                bpdict['bp_goterm'] = " ".join(bp[:-1])
                bpdict['bp_qoqid'] =  WDProp2QID_SPARQL(prop='P686', string= bpdict['bp_goid']).qid
                bpdict['bp_golabel'] = WDQID2Label_SPARQL(qid=bpdict['bp_qoqid']).label

                if bpdict['bp_goterm'] == bpdict['bp_golabel']:
                    bpgoqidlist.append(bpdict['bp_qoqid'])
                      #add log statement
            statements['biological_process'].append(PBB_Core.WDItemID(value=bpgoqidlist, prop_nr='P680', references=[copy.deepcopy(uniprot_protein_reference)]))


        final_statements = []

        for key in statements.keys():
            for statement in statements[key]:
                final_statements.append(statement)
                #print(key, statement.prop_nr, statement.value)

            wd_item_protein = PBB_Core.WDItemEngine(item_name=item_name, domain='proteins', data=final_statements)
            wd_item_protein.set_label(item_name)
            wd_item_protein.set_description(description)
            wd_item_protein.set_aliases(alias_list)
            pprint.pprint(wd_item_protein.get_wd_json_representation())
            #wd_item_protein.write(login)




login = PBB_login.WDLogin(sys.argv[1], sys.argv[2])
reference_genomes_list = NCBIReferenceGenomes()
for strain in reference_genomes_list.tid_list:
    mgi_record = MyGeneInfoRestBatchQuery(strain).gene_record
    unip_record = UniProtRESTBatchQuery(strain).enzyme_record


    combined = MGI_UNIP_Merger(mgi=mgi_record, unip=unip_record)

    for gene in combined.mgi_unip_dict:

        #wd_gene = WDGeneItem(gene, strain)
        #wd_gene.gene_item()
        wd_protein = WDProteinItem(gene, strain)
        wd_protein.protein_item()