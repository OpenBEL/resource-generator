#!/usr/bin/env python3

from operator import itemgetter
import datetime
import urllib.request
import os
import argparse
import annoheaders
from common import get_latest_MeSH_filename
from rdflib import Graph
from rdflib.namespace import RDF, RDFS, OWL, DC


def get_data(url):
    """ From url, get data, download and save locally. """
    REQ = urllib.request.urlopen(url)
    file_name = url.split('/')[-1]
    with open(file_name, 'b+w') as f:
        f.write(REQ.read())
    return file_name


def write_belanno(name, annodef, citation, anno_dict):
    """ Write .belanno file. """
    delim = '|'
    if len(anno_dict) == 0:
        print(
            '\tWARNING! No data found for {0}. Skipping writing .belanno file'.format(name))
    else:
        with open(name + '.belanno', 'w', encoding='utf-8') as belanno:
            belanno.write(annodef)
            belanno.write(annoheaders.author.substitute(year=year))
            belanno.write(citation + '\n')
            belanno.write(annoheaders.processing)
            belanno.write('[Values]\n')
            for k, v in sorted(anno_dict.items(), key=itemgetter(1)):
                belanno.write(v + delim + k + '\n')


def parse_owl(url, id, annotype):
    """" Download and parse owl files to return a dictionary of annotation values,
    where the unique identifier is the key and the name is the value,
    and the version and date information for the owl file. """
    # change working directory to 'source-data'
    if not os.path.exists('source-data'):
        os.mkdir('source-data')
    os.chdir('source-data')
    anno_dict = {}
    version = None
    owl = Graph()
    owl.parse(get_data(url))
    try:
        ontology = [s for s in owl.subjects(RDF.type, OWL.Ontology)][0]
    except:
        print("No Ontology URI for {0}".format(url))
    ver = owl.value(ontology, OWL.versionIRI)
    if ver:  # 1st try getting version and date from versionIRI
        version = str(ver)
        pub_date = version.split('/')[-2]
    if not version:  # 2nd try getting version from versionInfo
        ver = owl.value(ontology, OWL.versionInfo)
        if ver:
            version = str(ver[0])
        pub_date = owl.value(ontology, DC.date)
    if not pub_date:  # last try getting date from Ontology comments
        for o in owl.objects(ontology, RDFS.comment):
            if str(o).startswith('Date:'):
                pub_date = str(o).split(':')[-1].strip()
                pub_date = datetime.datetime.strptime(pub_date, '%dth %B %Y')
                pub_date = pub_date.strftime('%Y-%m-%d')

    print(version)
    print(pub_date)
    if id == 'EFO':  # only using cell_lines from EFO, make list of Classes
        cell_lines = [q[0] for q in owl.query(
            'SELECT * WHERE {?s rdfs:subClassOf+ <http://www.ebi.ac.uk/efo/EFO_0000322>}')]
    for s in owl.subjects(RDF.type, OWL.Class):
        val = owl.label(s)
        term = str(s).split('/')[-1]
        obsolete = owl.value(s, OWL.deprecated)
        if val and term.startswith(id) and not obsolete:
            if id == 'EFO' and s not in cell_lines:
                continue
            else:
                anno_dict[term] = val
    os.chdir(os.pardir)
    return anno_dict, ver, pub_date


def parse_mesh(mesh_url):
    """ Download and parse MeSH file, and return dictionary with
    UI (Unique identifier) as key and MH (MeSH Heading) and MNs (tree numbers)
    for each UI. """
    if not os.path.exists('source-data'):
        os.mkdir('source-data')
    os.chdir('source-data')
    with open(get_data(mesh_url), 'rt', encoding='utf-8') as mesh:
        UI = ""
        MH = ""
        MNs = set()
        MESH_dict = {}
        for line in iter(mesh):
            if not line.strip():
                continue
            if line.startswith('MH = '):
                MH = line.split('=')[1].strip()
            elif line.startswith('UI = '):
                UI = line.split('=')[1].strip()
            elif line.startswith('MN = '):
                MN = line.split('=')[1].strip()
                MNs.add(MN)
            elif line.startswith('*NEWRECORD'):
                # add UI, MH, MNs to dictionaries and set to empty
                if len(MNs) > 0:
                    MESH_dict[UI] = {'MH': MH, 'MNs': MNs}
                UI = ""
                MH = ""
                MNs = set()
    os.chdir(os.pardir)
    return MESH_dict

# data is dictionary key = name, value = list of tuples with source url and id
owl_data = {
    'anatomy': [('http://purl.obolibrary.org/obo/uberon.owl', 'UBERON')],
    'disease': [('http://purl.obolibrary.org/obo/doid.owl', 'DOID')],
    'cell': [('http://purl.obolibrary.org/obo/cl.owl', 'CL')],
    'cell-line': [('http://purl.obolibrary.org/obo/clo.owl', 'CLO'),
                  ('http://www.ebi.ac.uk/efo/efo.owl', 'EFO')]
}

# MeSH file is updated yearly - get latest
mesh_url = get_latest_MeSH_filename(
    'ftp://nlmpubs.nlm.nih.gov/online/mesh/.asciimesh/', 'd', '.bin')
# parse version from MeSH filename
mesh_ver = mesh_url.split('/')[-1].lstrip('d').rstrip('.bin')

# names of .belanno files from MeSH
mesh_anno_names = ['cell-structure', 'mesh-diseases', 'mesh-anatomy']

# MeSH sub-branches to imclude in mesh-anatomy
anatomy_branches = ('A01', 'A02', 'A03', 'A04', 'A05', 'A06', 'A07', 'A08',
                    'A09', 'A10', 'A11', 'A12', 'A14', 'A15', 'A16', 'A17')

disease_branches = ('C', 'F')

# command line argument
parser = argparse.ArgumentParser(
    description="""Generate BEL annotation files. """)
parser.add_argument(
    "-n",
    required=True,
    nargs=1,
    metavar="DIRECTORY",
    help="directory for new .belanno files")
args = parser.parse_args()

# date information for .belanno
today = datetime.datetime.today()
crdate = today.strftime('%Y-%m-%dT%T%Z')
version = today.strftime('%Y%m%d')
year = today.year

# create resoure directory and set to working dir
resource_dir = args.n[0]
if not os.path.exists(resource_dir):
    os.mkdir(resource_dir)
os.chdir(resource_dir)

# parse data from owl files and write .belanno files
for (name, dlist) in owl_data.items():
    anno_dict = {}
    citation = '[Citation]\n'
    print('\nGenerating .belanno file for {0} ...'.format(name))
    for (url, id) in dlist:
        (anno_values, sv, sdate) = parse_owl(url, id, name)
        anno_dict.update(anno_values)
        citation = citation + annoheaders.citation_info(id, sv, sdate)
    annodef = annoheaders.annotation_definition(name, version, crdate)
    write_belanno(name, annodef, citation, anno_dict)

# parse MeSH information and write MeSH .belanno files
MESH_dict = parse_mesh(mesh_url)
for anno in mesh_anno_names:
    anno_dict = {}
    print('Generating .belanno file for {0} ...'.format(anno))
    if anno == 'cell-structure':
        for UI, UI_dict in MESH_dict.items():
            if any('A11.284' in branch for branch in UI_dict['MNs']):
                anno_dict[UI] = UI_dict['MH']

    elif anno == 'mesh-diseases':
        for UI, UI_dict in MESH_dict.items():
            if any(branch.startswith(disease_branches)
                   for branch in UI_dict['MNs']):
                anno_dict[UI] = UI_dict['MH']

    elif anno == 'mesh-anatomy':
        for UI, UI_dict in MESH_dict.items():
            if any(branch.startswith(anatomy_branches)
                   for branch in UI_dict['MNs']):
                anno_dict[UI] = UI_dict['MH']
    annodef = annoheaders.annotation_definition(anno, version, crdate)
    citation = '[Citation]\n' + annoheaders.citation_info('MESH', mesh_ver)
    write_belanno(anno, annodef, citation, anno_dict)

