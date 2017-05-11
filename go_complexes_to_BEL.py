#!/usr/bin/env python3.3
from collections import defaultdict
import urllib.request
import csv
import gzip
from string import Template
import datetime
from bel_functions import bel_term

'''
Generates resource 'scaffolding' document for named complexes to be used with the BEL Compiler.

1. Identifies GOCC ID values that have the encoding 'C' (complex).
2. Downloads current GO Annotation files for human, mouse, and rat.
3. Validates values in GO annotation files as proteins ('GRP' encoding in .belns files).
4. Generates output file with complex hasComponent BEL statements.

'''

# name of output BEL file
output_file = 'named-complexes.bel'

# url for .belns document with complex, values are expected to be GO IDs
gocc_url = 'http://resources.openbel.org/belframework/testing/namespace/go-cellular-component-ids.belns'

# GO annotation file locations and associated namespaces
datasets = [
    ('9606',
     'ftp://ftp.geneontology.org/pub/go/gene-associations/gene_association.goa_human.gz',
     'HGNC'),
    ('10116',
     'ftp://ftp.geneontology.org/pub/go/gene-associations/gene_association.mgi.gz',
     'MGI'),
    ('10090',
     'ftp://ftp.geneontology.org/pub/go/gene-associations/gene_association.rgd.gz',
     'RGD')]


def get_complexes(belns_url):
    """ Parse .belns document and make set of complexes
    based on encoding values. Returns set of values that
    can be complexes and the keyword (namespace prefix)
    from the .belns document. """
    complexes = set()
    complex_ns = ''
    field = ''
    belns = urllib.request.urlopen(belns_url)
    for line in belns:
        line = line.decode('ISO-8859-1')
        if not line.strip():
            continue
        if line.startswith('['):
            field = line.strip()
            continue
        elif '[Namespace]' in field:
            if line.startswith('Keyword='):
                complex_ns = line.split('=')[1].strip()
        elif '[Values]' in field:
            (value, encoding) = line.split('|')
            encoding = encoding.strip()
            if encoding == 'C':
                complexes.add(value)
            else:
                continue
    return complexes, complex_ns


def get_encoding_dict(belns_url):
    """ Parse .belns document and return a value:encoding
    dictionary. Used to validate values as proteins. """
    field = ''
    encoding_dict = {}
    belns = urllib.request.urlopen(belns_url)
    for line in belns:
        line = line.decode('ISO-8859-1')
        if not line.strip():
            continue
        if line.startswith('['):
            field = line.strip()
            continue
        elif '[Values]' in field:
            (value, encoding) = line.split('|')
            encoding = encoding.strip()
            encoding_dict[value] = encoding
    return encoding_dict

# get complex ids from .belns file
complexes, complex_ns = get_complexes(gocc_url)
print('\nFound {0} complexes in {1}.'.format(len(complexes), gocc_url))

# info for bel script output file generation
annotations = {
    'Species': 'http://resources.openbel.org/belframework/testing/annotation/species-taxonomy-id.belanno'}
namespaces = {
    'HGNC': 'http://resources.openbel.org/belframework/testing/namespace/hgnc-approved-symbols.belns',
    'MGI': 'http://resources.openbel.org/belframework/testing/namespace/mgi-approved-symbols.belns',
    'RGD': 'http://resources.openbel.org/belframework/testing/namespace/rgd-approved-symbols.belns',
    complex_ns: gocc_url}

today = datetime.date.today()
version = str(today.year) + str(today.month) + str(today.day)
separator = '#' * 50

# generate file header
with open(output_file, 'w') as bel:
    bel.write(separator)
    bel.write('\n# Document Properties Section\n')
    bel.write('SET DOCUMENT Name = "GOCC Named Complex Scaffolding"\n')
    bel.write('SET DOCUMENT Description = "Named Complex scaffolding for use with the BEL compiler. File is generated from GO gene association files for human, mouse and rat."\n')
    bel.write('SET DOCUMENT Version = "{0}"\n'.format(version))
    bel.write('SET DOCUMENT Copyright = "Copyright (c) {0}, OpenBEL Project. This work is licensed under a Creative Commons Attribution 3.0 Unported License."\n'.format(
        str(today.year)))
    bel.write('SET DOCUMENT Authors = "OpenBEL"\n')
    bel.write('\n' + separator + '\n')
    bel.write('# Definitions Section\n')
    for ns, url in namespaces.items():
        bel.write('DEFINE NAMESPACE {0} AS URL "{1}"\n'.format(ns, url))
    for anno, url in annotations.items():
        bel.write('DEFINE ANNOTATION {0} AS URL "{1}"\n'.format(anno, url))
    bel.write(separator + '\n')
    bel.write('# Statements Section\n\n')

# download and open human, mouse, and rat gene association files from GO (datasets)
# write BEL file for each species with genes annotated to each complex term
for (species, url, ns) in datasets:
    # build dict to validate symbols as proteins
    encoding_dict = get_encoding_dict(namespaces.get(ns))
    date = ''
    # TODO - save to output/datasets directory
    print('Downloading {0} GO gene-associations data ...'.format(ns))
    go = urllib.request.urlopen(url)
    goa_file_name = url.split('/')[-1]
    with open(goa_file_name, 'b+w') as f:
        f.write(go.read())
    # make dictionary genemap with list of proteins for each complex GO term
    gomap = defaultdict(set)
    with gzip.open(goa_file_name, mode='rt') as g:
        reader = csv.reader(g, delimiter='\t')
        for row in reader:
            if row[0].startswith('!'):
                if 'Submission Date:' in row:
                    date = row[0].split(':')[1].strip()
                    (month, day, year) = date.split('/')
                    date = Template(
                        '${year}-${month}-${day}').substitute(year=year, day=day, month=month)
                else:
                    continue
            # find records with CC GO terms, that do not have a qualifier
            else:
                goid = row[4].replace('GO:', '')
                qualifier = row[3]
                symbol = row[2]
                if qualifier == '' and goid in complexes:
                    # convert complex into BEL term
                    term = bel_term(
                        row[4].replace(
                            'GO:',
                            ''),
                        complex_ns,
                        'complex')
                    # convert symbol into BEL term, exclude non-protein symbols
                    if encoding_dict.get(symbol) == 'GRP':
                        gene = bel_term(symbol, ns, 'p')
                        gomap[term].add(gene)

    # sort genes annotated to each complex in gomap
    gomap = {k: sorted(list(v)) for k, v in gomap.items()}
    print(
        '\tGenerating {0} BEL hasComponents statements.\n'.format(
            len(gomap)))
    # create statements and associated annotations; append to output file
    with open(output_file, 'a') as bel:
        bel.write('SET Species = {0}\n'.format(species))
        bel.write(
            'SET Citation = {3}"Online Resource", "GO Annotation File - {0}", "{1}", "{2}", "", ""{4}\n\n'.format(
                species, url, date, '{', '}'))
        for k, v in gomap.items():
            bel.write('{0} hasComponents list({1})\n'.format(k, ','.join(v)))
            # bel.write(k + ' hasComponents list(' + ",".join(v) + ')')
            # bel.write('\n')

