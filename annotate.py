# coding utf-8
#

'''
 annotate.py

 Python module to build three .belanno files.
   1. mesh-cell-structure.belanno
   2. mesh-diseases.belanno
   3. mesh-anatomy.belanno

'''

already_seen = set()


def make_annotations(d):
    delim = '|'
    branches = ['A01', 'A02', 'A03', 'A04', 'A05', 'A06', 'A07', 'A08',
                'A09', 'A10', 'A12', 'A14', 'A15', 'A16', 'A17']

    with open('mesh-cell-structure.belanno', 'w') as meshc, \
            open('mesh-diseases.belanno', 'w') as meshd, \
            open('mesh-anatomy.belanno', 'w') as mesha:
        for vals in d.get_annot_values():
            ui, mh, mns, synonyms = vals

            if any('A11.284' in branch for branch in mns):
                meshc.write(delim.join((mh, ui)) + '\n')

            elif any('C' in branch for branch in mns):
                meshd.write(delim.join((mh, ui)) + '\n')

            else:
                for branch in mns:
                    if any(b in branch for b in branches):
                        mesha.write(delim.join((mh, ui)) + '\n')
                        break

