#!/usr/bin/env python3.3
from string import Template
import datetime

namespaces = {'HGNC': 'hgnc-human-genes.belns',
              'MGI': 'mgi-mouse-genes.belns',
              'RGD': 'rgd-rat-genes.belns'}
annotations = {}
base_url = 'http://resources.openbel.org/belframework/latest/'
today = datetime.date.today()
version = today.strftime('%Y%m%d')


def bel_term(value, ns, f):
    """ Create bel term given value, namespace id,
    and bel function string. """
    must_quote_values = ['a', 'SET']
    must_quote_chars = [':', '(', ')', '<', '>', '.', '-', '/', '@', ' ']
    if any(char in value for char in must_quote_chars) or value in must_quote_values:
        s = Template('${f}(${ns}:"${value}")')
    else:
        s = Template('${f}(${ns}:${value})')
    term = s.substitute(f=f, ns=ns, value=value)
    return term


def write_bel_header(
        f,
        *,
        authors='OpenBEL',
        contact_info=None,
        doc_name="Document Name",
        description=None,
        licenses=None,
        version=version,
        namespaces=namespaces,
        annotations=annotations,
        base_url=base_url):
    '''Write BEL document header to file object f, given a document name,
     description, namespace dictionary (prefixes to urls) and annotation
    dictionary.'''
    separator = 50 * '#' + '\n'
    f.write(separator + '# Document Properties Section\n')
    f.write('SET DOCUMENT Name = "{0}"\n'.format(doc_name))
    if description:
        f.write('SET DOCUMENT Description = "{0}"\n'.format(description))
    f.write('SET DOCUMENT Version = "{0}"\n'.format(version))
    f.write('SET DOCUMENT Copyright = "Copyright (c) {0}, OpenBEL Project. This work is licensed under a Creative Commons Attribution 3.0 Unported License."\n'.format(
        str(today.year)))
    f.write('SET DOCUMENT Authors = {0}\n'.format(authors))
    if licenses:
        f.write('SET DOCUMENT Licenses = {0}\n'.format(licenses))
    if contact_info:
        f.write('SET Document ContactInfo = {0}\n'.format(contact_info))
    f.write('\n' + separator + '# Definitions Section\n')
    for ns_prefix, ns_name in namespaces.items():
        base_url = base_url.rstrip('/')
        f.write(
            'DEFINE NAMESPACE {0} AS URL "{1}/namespace/{2}"\n'.format(ns_prefix, base_url, ns_name))
    f.write('\n')
    for anno, url in annotations.items():
        f.write(
            'DEFINE ANNOTATION {0} AS URL "{1}/annotation/{2}"\n'.format(anno, base_url, url))
    f.write(separator)
    f.write('# Statements Section\n\n')
    return None

