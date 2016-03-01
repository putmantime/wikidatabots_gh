__author__ = 'andra'

#!/usr/bin/python
"""
This example shows basic document generation functionality.
..  :copyright: (c) 2014 by Jelte Fennema.
    :license: MIT, see License for more details.
"""

# begin-doc-include
from pylatex import Document, Section, Subsection, Command, Tabular, Figure, Package, TikZ, Axis, Plot
from pylatex.utils import italic, NoEscape
import pandas as pd
import matplotlib.pyplot as plt  # noqa


if __name__ == '__main__':
    # Basic document
    doc = Document('basic')

    doc.generate_pdf()
    doc.generate_tex()

    # Document with `\maketitle` command activated
    doc = Document()
    doc.packages.append(Package('geometry', options=['left=2cm', 'right=2cm']))

    doc.preamble.append(Command('title', 'Curation report on Diseasebot'))
    doc.preamble.append(Command('author', 'Andra Waagmeester'))
    doc.preamble.append(Command('date', NoEscape(r'\today')))
    doc.append(NoEscape(r'\maketitle'))

    doc.generate_pdf('basic_maketitle', clean=False)

    df = pd.read_csv('/tmp/logs/WD_bot_run-2016-02-24_16:10.log', sep=',', quotechar='"', encoding='latin-1', header=None, error_bad_lines=False)
    print(df[0].value_counts())

    with doc.create(Section('Items processed')):
        counts = df[0].value_counts()
        table1 = Tabular('|c|c|')
        table1.add_hline()
        table1.add_row(("INFO", "ERROR"))
        table1.add_row((counts.INFO, counts.ERROR))
        table1.add_hline()
        doc.append(table1)

    with doc.create(Section('Error types and counts')):
        error_df = df[df[0] == "ERROR"]
        errorType_counts = error_df[3].value_counts()
        errorType_details = error_df[4].value_counts()
        error_table = Tabular('|c|c|')
        error_table.add_hline()
        error_table.add_row(("ERROR_TYPE", "counts"))
        error_table.add_hline()

    with doc.create(Section('I am a section')):
        doc.append('Take a look at this beautiful plot:')

        with doc.create(Figure(position='htbp')) as plot:
            plot.add_plot(width=NoEscape('200pt'))
            plot.add_caption('I am a caption.')

        doc.append('Created using matplotlib.')



    # Add stuff to the document
    table2 = Tabular('|c|c|c|c|c|c|c|')
    table2.add_hline()


    with doc.create(Section('Appendix X: Full log file')):
        for row in df.itertuples():
            table2.add_row((row[1], row[2], row[3], row[4], row[5], row[6], row[7]))
            table2.add_hline()
            #doc.append(row[1])
            #print(row)
    doc.append(table2)
    doc.generate_pdf('/tmp/basic_maketitle2')
    tex = doc.dumps()  # The document as string in LaTeX syntax