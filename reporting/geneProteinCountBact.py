
__author__ = 'tim'

import sys
print(sys.version)
from SPARQLWrapper import SPARQLWrapper, JSON
import pprint
from pylatex import Document, Section, Subsection, Command, Tabular, Figure, Package, TikZ, Axis, Plot
from pylatex.utils import italic, NoEscape

def runGeneCounts(countDict):
    sparql = SPARQLWrapper("https://query.wikidata.org/bigdata/namespace/wdq/sparql")
    sparql.setQuery("""
        PREFIX wd: <http://www.wikidata.org/entity/>
        PREFIX wdt: <http://www.wikidata.org/prop/direct/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

        SELECT ?species ?label (count(distinct ?entrezID) as ?gene_counts) WHERE {
           ?gene wdt:P351 ?entrezID ; # P351 Entrez Gene ID
                 wdt:P703 ?species . # P703 Found in taxon
           # ?protein wdt:P352 ?uniprotID ;
           #         wdt:P703 ?species .
           ?species wdt:P171* wd:Q10876 ;
                rdfs:label ?label .
           filter (lang(?label) = "en") .
         }
         GROUP BY ?species ?label
        """)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    print(len(results["results"]["bindings"]))
    pprint.pprint(results["results"]["bindings"])
    for result in results["results"]["bindings"]:
        if result["label"]["value"] not in countDict.keys():
            countDict[result["label"]["value"]] = dict()
        countDict[result["label"]["value"]]["gene_counts"] = result["gene_counts"]["value"]
        countDict[result["label"]["value"]]["wd_uri"] = result["species"]["value"]

def runProteinCounts(countDict):
    sparql = SPARQLWrapper("https://query.wikidata.org/bigdata/namespace/wdq/sparql")
    sparql.setQuery("""
        PREFIX wd: <http://www.wikidata.org/entity/>
        PREFIX wdt: <http://www.wikidata.org/prop/direct/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

        SELECT ?species ?label (count(distinct ?uniprotID) as ?protein_counts) WHERE {
           #?gene wdt:P351 ?entrezID ; # P351 Entrez Gene ID
           #      wdt:P703 ?species . # P703 Found in taxon
            ?protein wdt:P352 ?uniprotID ;
                    wdt:P703 ?species .
           ?species wdt:P171* wd:Q10876 ;
                rdfs:label ?label .
           filter (lang(?label) = "en") .
         }
         GROUP BY ?species ?label
        """)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    print(len(results["results"]["bindings"]))
    pprint.pprint(results["results"]["bindings"])
    for result in results["results"]["bindings"]:
        if result["label"]["value"] not in countDict.keys():
            countDict[result["label"]["value"]] = dict()
        countDict[result["label"]["value"]]["protein_counts"] = result["protein_counts"]["value"]
        countDict[result["label"]["value"]]["wd_uri"] = result["species"]["value"]


counts = dict()
runGeneCounts(counts)
runProteinCounts(counts)
pprint.pprint(counts)


# Basic document
doc = Document('basic')
doc.packages.append(Package('hyperref'))
doc.preamble.append(Command('title', 'Number of genes and proteins per bacterial strain'))
doc.preamble.append(Command('author', 'Tim Putman'))
doc.preamble.append(Command('date', NoEscape(r'\today')))
doc.append(NoEscape(r'\maketitle'))
doc.generate_pdf('basic_maketitle', clean=False)

with doc.create(Section('Genes')):
    table1 = Tabular('|l|l|l|}')
    table1.add_hline()
    table1.add_row(("species", "genes", "proteins"))
    table1.add_hline()

    for species in counts.keys():
        if "gene_counts" not in counts[species].keys():
            table1.add_row((Command('href', options=None,arguments=counts[species]["wd_uri"], extra_arguments=species), "0", counts[species]["protein_counts"]))
        elif "protein_counts" not in counts[species].keys():
            table1.add_row((Command('href', options=None,arguments=counts[species]["wd_uri"], extra_arguments=species), counts[species]["gene_counts"], "0))
        else:
            table1.add_row((Command('href', options=None,arguments=counts[species]["wd_uri"], extra_arguments=species), counts[species]["gene_counts"], counts[species]["protein_counts"]))
    table1.add_hline()
    doc.append(table1)

speciesList = []
for species in counts.keys():
    speciesList.append(species)
with doc.create(Subsection('Beautiful graphs')):
    with doc.create(TikZ()):
            plot_options = "title = Genes and proteins per species, symbolic y coords = {"+speciesList.join(",")+"}"
            with doc.create(Axis(options=plot_options)) as plot:
                print(14)

            #with doc.create(Axis(options=plot_options)) as plot:

doc.generate_pdf('/tmp/microbio')
tex = doc.dumps()
pprint.pprint(tex)



