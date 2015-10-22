# coding: utf-8

'''
 parsers.py

 Contains parser objects constructed for various datasets that
 will be used to build the .belns, .beleq, and .belanno files.

'''

from common import gzip_to_text
from lxml import etree
import os
import csv
import gzip
import urllib.request
import zipfile
import io

class Parser(object):
    def __init__(self, url):
        self.url = url

    def parse():
        pass


class EntrezGeneInfoParser(Parser):
    resourceLocation = """http://resource.belframework.org/belframework/1.0/
                          namespace/entrez-gene-ids-hmr.belns"""

    def __init__(self, url):
        super(EntrezGeneInfoParser, self).__init__(url)
        self.entrez_info = url

    def parse(self):

        # columns for an Entrez gene info dataset.
        entrez_info_headers = ['tax_id', 'GeneID', 'Symbol', 'LocusTag',
                               'Synonyms', 'dbXrefs', 'chromosome',
                               'map_location', 'description',
                               'type_of_gene',
                               'Symbol_from_nomenclature_authority',
                               'Full_name_from_nomenclature_authority',
                               'Nomenclature_status',
                               'Other_designations', 'Modification_date']

        # dictionary for base gene info
        info_csvr = csv.DictReader(gzip_to_text(self.entrez_info),
                                   delimiter='\t',
                                   fieldnames=entrez_info_headers)

        for row in info_csvr:
            if row['tax_id'] in ('9606', '10090', '10116'):
                yield row

    def __str__(self):
        return "EntrezGeneInfo_Parser"


class EntrezGeneHistoryParser(Parser):
    resourceLocation = """"http://resource.belframework.org/belframework/1.0/
                           namespace/entrez-gene-ids-hmr.belns"""

    def __init__(self, url):
        super(EntrezGeneHistoryParser, self).__init__(url)
        self.entrez_history = url

    def parse(self):

        entrez_history_headers = ["tax_id", "GeneID", "Discontinued_GeneID",
                                  "Discontinued_Symbol", "Discontinue_Date"]

        # dictionary for base gene info
        history_csvr = csv.DictReader(gzip_to_text(self.entrez_history),
                                      delimiter='\t',
                                      fieldnames=entrez_history_headers)

        for row in history_csvr:
            if row['tax_id'] in ("9606", "10090", "10116"):
                yield row

    def __str__(self):
        return "EntrezGeneHistory_Parser"


class HGNCParser(Parser):
    resourceLocation = """http://resource.belframework.org/belframework/1.0/
                          namespace/hgnc-approved-symbols.belns"""

    def __init__(self, url):
        super(HGNCParser, self).__init__(url)
        self.hgnc_file = url

    def parse(self):

        # use iso-8859-1 as default encoding.
        with open(self.hgnc_file, "r", encoding="iso-8859-1") as hgncf:

            # Note that HGNC uses TWO columns named the same thing for Entrez
            # Gene ID. Currently we are not using these columns and it is not a
            # big deal, but in the future we could account for this by using
            # custom headers (like EntrezGeneInfo_Parser), or resolving to the
            # SECOND of the Entrez Gene ID columns.
            hgnc_csvr = csv.DictReader(hgncf, delimiter='\t')

            for row in hgnc_csvr:
                yield row

    def __str__(self):
        return "HGNC_Parser"


class MGIParser(Parser):
    resourceLocation = """http://resource.belframework.org/belframework/1.0/
                          namespace/mgi-approved-symbols.belns"""

    def __init__(self, url):
        super(MGIParser, self).__init__(url)
        self.mgi_file = url

    def parse(self):
        with open(self.mgi_file, "r") as mgif:
            mgi_csvr = csv.DictReader(mgif, delimiter='\t')

            for row in mgi_csvr:
                yield row

    def __str__(self):
        return "MGI_Parser"


class RGDParser(Parser):
    resourceLocation = """http://resource.belframework.org/belframework/1.0/
                          namespace/rgd-approved-symbols.belns"""

    def __init__(self, url):
        super(RGDParser, self).__init__(url)
        self.rgd_file = url

    def parse(self):
        with open(self.rgd_file, "r") as rgdf:
            # skip all the comment lines beginning with '#' and also the header.
            rgd_csvr = csv.DictReader(filter(lambda row:
                                                 not row[0].startswith('#'), rgdf),
                                      delimiter='\t')

            for row in rgd_csvr:
                yield row

    def __str__(self):
        return "RGD_Parser"


# This class exists mainly as a way to break the iteration loop during parsing
# of the SwissProt dataset if needed.
class GeneTypeError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class SwissProtParser(Parser):
    resourceLocation_accession_numbers = """http://resource.belframework.org/
                belframework/1.0/namespace/swissprot-accession-numbers.belns"""
    resourceLocation_entry_names = """http://resource.belframework.org/
                      belframework/1.0/namespace/swissprot-entry-names.belns"""

    def __init__(self, url):
        super(SwissProtParser, self).__init__(url)
        self.sprot_file = url
        self.entries = {}
        self.accession_numbers = {}
        self.gene_ids = {}
        self.tax_ids = {'9606', '10090', '10116'}
        self.pro = '{http://uniprot.org/uniprot}protein'
        self.rec_name = '{http://uniprot.org/uniprot}recommendedName'
        self.full_name = '{http://uniprot.org/uniprot}fullName'
        self.short_name = '{http://uniprot.org/uniprot}shortName'
        self.alt_name = '{http://uniprot.org/uniprot}alternativeName'
        self.db_ref = '{http://uniprot.org/uniprot}dbReference'
        self.organism = '{http://uniprot.org/uniprot}organism'
        self.entry = '{http://uniprot.org/uniprot}entry'
        self.accession = '{http://uniprot.org/uniprot}accession'
        self.name = '{http://uniprot.org/uniprot}name'

    def parse(self):

        with gzip.open(self.sprot_file) as sprotf:
#        with open(self.sprot_file, 'rb') as sprotf:
            ctx = etree.iterparse(sprotf, tag=self.entry)

            for ev, e in ctx:
                temp_dict = {}
                n_dict = {}

                # stop evaluating if this entry is not in the Swiss-Prot dataset
                if e.get('dataset') != 'Swiss-Prot':
                    e.clear()
                    continue

                # stop evaluating if this entry is not for human, mouse, or rat
                org = e.find(self.organism)

                # use a custom exception to break to next iteration (e)
                # if tax ref is not found.
                try:
                    for org_child in org:
                        if org_child.tag == self.db_ref:
                            # restrict by NCBI Taxonomy reference
                            if org_child.get('id') not in self.tax_ids:
                                e.clear()
                                raise GeneTypeError(org_child.get('id'))
                            else:
                                # add NCBI Taxonomy and the id for the entry
                                # to the dict
                                temp_dict[org_child.get('type')] = \
                                    org_child.get('id')
                except GeneTypeError:
                    continue

                # get entry name, add it to the dict
                entry_name = e.find(self.name).text
                temp_dict['name'] = entry_name

                # get protein data, add recommended full and short names to dict
                protein = e.find(self.pro)

                for child in protein.find(self.rec_name):
                    if child.tag == self.full_name:
                        temp_dict['recommendedFullName'] = child.text
                    if child.tag == self.short_name:
                        temp_dict['recommendedShortName'] = child.text
                alt_shortnames = []
                alt_fullnames = []

                protein = e.find(self.pro)
                for altName in protein.findall(self.alt_name):
                    for child in altName:
                        if child.tag == self.full_name:
                            alt_fullnames.append(child.text)
                        if child.tag == self.short_name:
                            alt_shortnames.append(child.text)

                temp_dict['alternativeFullNames'] = alt_fullnames
                temp_dict['alternativeShortNames'] = alt_shortnames

                # get all accessions
                entry_accessions = []
                for entry_accession in e.findall(self.accession):
                    acc = entry_accession.text
                    entry_accessions.append(acc)
                    if acc in self.accession_numbers:
                        self.accession_numbers[acc] = None
                    else:
                        self.accession_numbers[acc] = 1

                # add the array of accessions to the dict
                temp_dict["accessions"] = entry_accessions

                # add dbReference type (human, rat, and mouse) and gene ids to
                # the dict
                type_set = ['GeneId', 'MGI', 'HGNC', 'RGD']
                for dbr in e.findall(self.db_ref):
                    if dbr.get('type') in type_set:
                        gene_id = dbr.get('id')
                        if dbr.get('type') not in n_dict:
                            n_dict[dbr.get('type')] = [gene_id]
                        else:
                            n_dict[dbr.get('type')].append(gene_id)
                temp_dict['dbReference'] = n_dict

                # clear the tree before next iteration
                e.clear()
                while e.getprevious() is not None:
                    del e.getparent()[0]

                yield temp_dict

    def __str__(self):
        return 'SwissProt_Parser'


# Helper function for AffyParser. This will save each of the downloaded
# URLs and return the file pointer.
def get_data(url):
    # from url, download and save file
    REQ = urllib.request.urlopen(url)
    file_name = url.split('/')[-1]
    os.chdir('datasets/')
    with open(file_name,'wb') as f:
        f.write(REQ.read())
    os.chdir('../')
    return file_name

def filter_plus_print(row):
    return not row.startswith('#')


class AffyParser(Parser):

    def __init__(self, url):
        super(AffyParser, self).__init__(url)
        self.affy_file = url

    def parse(self):

        # the arrays we are concerned with
        array_names = ['HG-U133A', 'HG-U133B', 'HG-U133_Plus_2', 'HG_U95Av2',
                       'MG_U74A', 'MG_U74B', 'MG_U74C', 'MOE430A', 'MOE430B',
                       'Mouse430A_2', 'Mouse430_2', 'RAE230A', 'RAE230B',
                       'Rat230_2']

        urls = []
        with open(self.affy_file, 'rb') as affyf:
            ctx = etree.iterparse(affyf, events=('start', 'end'))

            # This is probably not the best way to traverse this tree. Look at
            # the lxml.etree API more closely for possible implementations when
            # refactoring
            # NOTES - put some debugging in here to see how this is parsing,
            # may be a better way to parse (like using diff events).
            for ev, e in ctx:
                # iterate the Array elements
                for n in e.findall('Array'):
                    name = n.get('name')
                    if name in array_names:
                        # iterate Annotation elements
                        for child in n:
                            if child.get('type') == 'Annot CSV':
                                # iterate File elements
                                for g_child in child:
                                    # get the URL and add to the list
                                    for gg_child in g_child:
                                        urls.append(gg_child.text)

        # iterate over the list of URLs returned from the Affy XML feed
        for link in urls:
            affy_reader = {}

            # get_data() downloads the file, saves it as a .csv.zip, and
            # returns a pointer to the file.
            n = get_data(link)
            z = zipfile.ZipFile('datasets/'+n, 'r')

            # only want the .csv from the archive (also contains a .txt)
            for name in z.namelist():
                if '.csv' in name:
                    print('\tExtracting - ' +name)
                    # wrap in a TextIOWrapper. otherwise it returns bytes.
                    affy_reader = csv.DictReader(filter(lambda x:
                                                        not x.startswith('#'),
                                                        io.TextIOWrapper(z.open(name))),
                                                        delimiter=',')

                    for x in affy_reader:
                        yield x

    def __str__(self):
        return 'Affy_Parser'

class Gene2AccParser(Parser):

    def __init__(self, url):
       super(Gene2AccParser, self).__init__(url)
       self.gene2acc_file = url

    def parse(self):

        # would like to have DictReader handle this, but need a way to
        # deal with the special case of the first value beginning with
        # a hashtag. i.e. #Format: <-- is NOT a column header.
        column_headers = ['tax_id', 'GeneID', 'status',
                          'RNA nucleotide accession.version',
                          'RNA nucleotide gi', 'protein accession.version',
                          'protein gi', 'genomic nucleotide accession.version',
                          'genomic nucleotide gi',
                          'start position on the genomic accession',
                          'end position on the genomic accession',
                          'orientation', 'assembly',
                          'mature peptide accession.version',
                          'mature peptide gi', 'Symbol']

        g2a_reader = csv.DictReader(gzip_to_text(self.gene2acc_file), delimiter='\t',
                                    fieldnames=column_headers)

        for row in g2a_reader:
            yield row

    def __str__(self):
        return 'Gene2Acc_Parser'

class BELNamespaceParser(Parser):

    def __init__(self):
        self.old_files = 'http://resource.belframework.org./belframework/1.0/index.xml'
        self.anno_def = '{http://www.belscript.org/schema/annotationdefinitions}annotationdefinitions'
        self.namespace = '{http://www.belscript.org/schema/namespace}namespace'
        self.namespaces = '{http://www.belscript.org/schema/namespaces}namespaces'

    def parse(self):

        tree = etree.parse(self.old_files)

        # xpath will return all elements under this namespace (list of bel namespace urls)
        urls = tree.xpath('//*[local-name()="namespace"]/@idx:resourceLocation',
                          namespaces={'xsi': 'http://www.w3.org/2001/XMLSchema-instance',
                                      'idx' : 'http://www.belscript.org/schema/index'})

        for url in urls:
            yield url

    def __str__(self):
        return 'BELNamespace_Parser'


class BELEquivalenceParser(Parser):

    def __init__(self):
        self.old_files = 'http://resource.belframework.org./belframework/1.0/index.xml'
        self.anno_def = '{http://www.belscript.org/schema/annotationdefinitions}annotationdefinitions'
        self.namespace = '{http://www.belscript.org/schema/namespace}namespace'
        self.namespaces = '{http://www.belscript.org/schema/namespaces}namespaces'

    def parse(self):

        tree = etree.parse(self.old_files)

        # xpath will return all elements under this namespace (list of bel equivalence urls)
        urls = tree.xpath('//*[local-name()="equivalence"]/@idx:resourceLocation',
                          namespaces={'xsi': 'http://www.w3.org/2001/XMLSchema-instance',
                                      'idx' : 'http://www.belscript.org/schema/index'})

        for url in urls:
            yield url

    def __str__(self):
        return 'BELEquivalence_Parser'


class BELAnnotationsParser(Parser):

    def __init__(self):
        self.old_files = 'http://resource.belframework.org./belframework/1.0/index.xml'

    def parse(self):

        tree = etree.parse(self.old_files)

        # xpath will return all elements under this namespace (list of bel equivalence urls)
        urls = tree.xpath('//*[local-name()="annotationdefinition"]/@idx:resourceLocation',
                          namespaces={'xsi': 'http://www.w3.org/2001/XMLSchema-instance',
                                      'idx' : 'http://www.belscript.org/schema/index'})

        for url in urls:
            yield url

    def __str__(self):
        return 'BELAnnotations_Parser'

# This one uses iterparse(), much faster than xpath on the
# bigger .owl file.
class CHEBIParser(Parser):

    def __init__(self, url):
        super(CHEBIParser, self).__init__(url)
        self.chebi_file = url
        self.classy = '{http://www.w3.org/2002/07/owl#}Class'
        self.label = '{http://www.w3.org/2000/01/rdf-schema#}label'
        self.altId = '{http://purl.obolibrary.org/obo#}altId'
        self.synonym = '{http://purl.obolibrary.org/obo#}Synonym'

    def parse(self):

        with open(self.chebi_file, 'rb') as cf:
            tree = etree.iterparse(cf, tag=self.classy)
            for event, elem in tree:
                if len(elem.values()) != 0:
                    chebi_dict = {}
                    synonyms = set()
                    alt_ids = set()
                    name = ''
                    vals = elem.values()
                    chebi_dict['primary_id'] = vals[0].split('CHEBI_')[1]
                    children = elem.getchildren()
                    for child in children:
                        if child.tag == self.label:
                            name = child.text
                        if child.tag == self.altId:
                            alt_ids.add(child.text.split(':')[1])
                        if child.tag == self.synonym:
                            synonyms.add(child.text)
                        chebi_dict['name'] = name
                        chebi_dict['alt_ids'] = alt_ids
                        chebi_dict['synonyms'] = synonyms

                    yield chebi_dict

    def __str__(self):
        return 'CHEBI_Parser'


class PubNamespaceParser(Parser):

    def __init__(self, url):
       super(PubNamespaceParser, self).__init__(url)
       self.pub_file = url

    def parse(self):
        column_headers = ['pubchem_id', 'synonym']
        pub_reader = csv.DictReader(gzip_to_text(self.pub_file), delimiter='\t',
                                    fieldnames=column_headers)

        with open(self.pub_file, 'r') as fp:
            pub_reader = csv.DictReader(fp, delimiter='\t',
                                        fieldnames=column_headers)

        for row in pub_reader:
            yield row

    def __str__(self):
        return 'PubNamespace_Parser'


class PubEquivParser(Parser):

    def __init__(self, url):
       super(PubEquivParser, self).__init__(url)
       self.cid_file = url

    def parse(self):

        column_headers = ['PubChem SID', 'Source', 'External ID', 'PubChem CID']
        cid_reader = csv.DictReader(gzip_to_text(self.cid_file), delimiter='\t',
                                    fieldnames=column_headers)

        for row in cid_reader:
            yield row

    def __str__(self):
        return 'PubEquiv_Parser'


class SCHEMParser(Parser):

    def __init__(self, url):
        super(SCHEMParser, self).__init__(url)
        self.schem_file = url

    def parse(self):
        isFalse = True
        with open(self.schem_file, 'r') as fp:
            for line in fp.readlines():
                if '[Values]' not in line and isFalse:
                    continue
                elif '[Values]' in line:
                    isFalse = False
                    continue
                else:
                    schem_id = line.split('|')[0]
                    yield {'schem_id' : schem_id }

    def __str__(self):
        return 'SCHEM_Parser'


class SDISParser(Parser):

    def __init__(self, url):
        super(SDISParser, self).__init__(url)
        self.sdis_file = url

    def parse(self):
        isFalse = True
        with open(self.sdis_file, 'r') as fp:
            for line in fp.readlines():
                if '[Values]' not in line and isFalse:
                    continue
                elif '[Values]' in line:
                    isFalse = False
                    continue
                else:
                    sdis_id = line.split('|')[0]
                    yield {'sdis_id' : sdis_id }

    def __str__(self):
        return 'SDIS_Parser'


class SCHEMtoCHEBIParser(Parser):

    def __init__(self, url):
        super(SCHEMtoCHEBIParser, self).__init__(url)
        self.schem_to_chebi = url

    def parse(self):

        column_headers = ['SCHEM_term', 'CHEBIID', 'CHEBI_name']

        with open(self.schem_to_chebi, 'r') as fp:
            rdr = csv.DictReader(fp, delimiter='\t', fieldnames=column_headers)

            for row in rdr:
                yield row

    def __str__(self):
        return 'SCHEMtoCHEBI_Parser'


class SDIStoDOParser(Parser):

    def __init__(self, url):
        super(SDIStoDOParser, self).__init__(url)
        self.sdis_to_do = url

    def parse(self):

        column_headers = ['SDIS_term', 'DOID', 'DO_name']

        with open(self.sdis_to_do, 'r') as fp:
            rdr = csv.DictReader(fp, delimiter='\t', fieldnames=column_headers)

            for row in rdr:
                yield row

    def __str__(self):
        return 'SDIStoDO_Parser'


class GOBPParser(Parser):

    def __init__(self, url):
        super(GOBPParser, self).__init__(url)
        self.go_file = url

    def parse(self):

        # parse xml tree using lxml
        parser = etree.XMLParser(ns_clean=True, recover=True, encoding='UTF-8')
        with gzip.open(self.go_file, 'r') as go:
            root = etree.parse(go, parser)
            bp_terms = root.xpath("/obo/term [namespace = 'biological_process' and not(is_obsolete)]")

            # iterate the NON-OBSOLETE biological_process terms
            for t in bp_terms:
                bp_termid = t.find('id').text
                bp_termname = t.find('name').text
                if t.findall('alt_id') is not None:
                    bp_altids = [x.text for x in t.findall('alt_id')]
                else:
                    bp_altids = False
                yield { 'termid' : bp_termid, 'termname' : bp_termname,
                        'alt_ids' : bp_altids }

    def obsolete_parse(self):

        # parse xml tree using lxml
        parser = etree.XMLParser(ns_clean=True, recover=True, encoding='UTF-8')
        with gzip.open(self.go_file, 'r') as go:
            root = etree.parse(go, parser)
            bp_terms = root.xpath("/obo/term [namespace = 'biological_process' and is_obsolete]")

            # iterate the OBSOLETE biological_process terms
            for t in bp_terms:
                bp_termid = t.find('id').text
                bp_termname = t.find('name').text
                if t.findall('alt_id') is not None:
                    bp_altids = [x.text for x in t.findall('alt_id')]
                else:
                    bp_altids = False
                yield { 'termid' : bp_termid, 'termname' : bp_termname,
                        'alt_ids' : bp_altids }


    def __str__(self):
        return 'GOBP_Parser'


class GOCCParser(Parser):

    def __init__(self, url):
        super(GOCCParser, self).__init__(url)
        self.go_file = url
        self.mesh_file = 'meshcs_to_gocc.csv'

    def parse(self):

        # initialize empty dictionaries using tuple assignment
        cc_parents, accession_dict, term_dict = {}, {}, {}

        # parse xml tree using lxml
        parser = etree.XMLParser(ns_clean=True, recover=True, encoding='UTF-8')
        root = etree.parse(self.go_file, parser)
        cc_terms = root.xpath("/obo/term [namespace = 'cellular_component' and not(is_obsolete)]")

        # iterate the NON-OBSOLETE complex terms to build parent dictionary
        for t in cc_terms:
            cc_termid = t.find("id").text
            cc_parent_ids = [isa.text for isa in t.findall("is_a")]
            cc_parents[cc_termid] = cc_parent_ids

        for t in cc_terms:
            cc_termid = t.find('id').text
            cc_termname = t.find('name').text
            cc_parent_stack = cc_parents[cc_termid]
            if t.findall('alt_id') is not None:
                cc_altids = [x.text for x in t.findall('alt_id')]
            else:
                cc_altids = False
            complex = False

            if cc_termid == "GO:0032991":
                complex = True
            elif t.find("is_root") is not None:
                complex = False
            else:
                cc_parent_stack.extend(cc_parents[cc_termid])
                while len(cc_parent_stack) > 0:
                    cc_parent_id = cc_parent_stack.pop()

                    if cc_parent_id == "GO:0032991":
                        complex = True
                        break

                    if cc_parent_id in cc_parents:
                        cc_parent_stack.extend(cc_parents[cc_parent_id])

            yield { 'termid' : cc_termid, 'termname' : cc_termname,
                    'altids' : cc_altids, 'complex' : complex }

    def obsolete_parse(self):

        # initialize empty dictionaries using tuple assignment
        cc_parents, accession_dict, term_dict = {}, {}, {}

        # parse xml tree using lxml
        parser = etree.XMLParser(ns_clean=True, recover=True, encoding='UTF-8')
        root = etree.parse(self.go_file, parser)
        cc_terms = root.xpath("/obo/term [namespace = 'cellular_component' and is_obsolete]")

        # iterate the OBSOLETE complex terms to build parent dictionary
        for t in cc_terms:
            cc_termid = t.find("id").text
            cc_parent_ids = [isa.text for isa in t.findall("is_a")]
            cc_parents[cc_termid] = cc_parent_ids

        for t in cc_terms:
            cc_termid = t.find('id').text
            cc_termname = t.find('name').text
            cc_parent_stack = cc_parents[cc_termid]
            if t.findall('alt_id') is not None:
                cc_altids = [x.text for x in t.findall('alt_id')]
            else:
                cc_altids = False
            complex = False

            if cc_termid == "GO:0032991":
                complex = True
            elif t.find("is_root") is not None:
                complex = False
            else:
                cc_parent_stack.extend(cc_parents[cc_termid])
                while len(cc_parent_stack) > 0:
                    cc_parent_id = cc_parent_stack.pop()

                    if cc_parent_id == "GO:0032991":
                        complex = True
                        break

                    if cc_parent_id in cc_parents:
                        cc_parent_stack.extend(cc_parents[cc_parent_id])

            yield { 'termid' : cc_termid, 'termname' : cc_termname,
                    'altids' : cc_altids, 'complex' : complex }

    def __str__(self):
        return 'GOCC_Parser'


class MESHParser(Parser):

    def __init__(self, url):
        super(MESHParser, self).__init__(url)
        self.mesh_file = url

    def parse(self):

        # ui - unique identifier / mh - mesh header
        # mn - tree # / st - semantic type
        ui = ''
        mh = ''
        mns = set()
        sts = set()
        synonyms = set()
        firstTime = True
        with open(self.mesh_file, 'r') as fp:
            for line in fp.readlines():
                if line.startswith('MH ='):
                    mh = line.split('=')[1].strip()
                elif line.startswith('UI ='):
                    ui = line.split('=')[1].strip()
                elif line.startswith('MN ='):
                    mn = line.split('=')[1].strip()
                    mns.add(mn)
                elif line.startswith('ST ='):
                    st = line.split('=')[1].strip()
                    sts.add(st)
                elif line.startswith('PRINT ENTRY ='):
                    entry = line.split('=')[1].strip()
                    if '|EQV|' in line:
                        entries = entry.split('|')
                        last = entries[-1]
                        num_syns = last.count('a')
                        while num_syns > 0:
                            num_syns = num_syns - 1
                            s = entries[num_syns]
                            synonyms.add(s.strip())
                    else:
                        if '|' not in entry:
                            synonyms.add(entry)
                elif line.startswith('ENTRY ='):
                    entry = line.split('=')[1].strip()
                    if '|EQV|' in line:
                        entries = entry.split('|')
                        last = entries[-1]
                        num_syns = last.count('a')
                        while num_syns > 0:
                            num_syns = num_syns - 1
                            s = entries[num_syns]
                            synonyms.add(s.strip())
                    else:
                        if '|' not in entry:
                            synonyms.add(entry)
                elif line.startswith('*NEWRECORD'):
                    # file begins with *NEWRECORD so skip that one (dont yield)
                    if firstTime:
                        firstTime = False
                        continue
                    else:

                        yield { 'ui' : ui, 'mesh_header' : mh,
                                'mns' : mns, 'sts' : sts,
                                'synonyms' : synonyms }
                        ui = ''
                        mh = ''
                        mns = set()
                        sts = set()
                        synonyms = set()

    def __str__(self):
        return 'MESH_Parser'


class SwissWithdrawnParser(Parser):

    def __init__(self, url):
        super(SwissWithdrawnParser, self).__init__(url)
        self.s_file = url

    def parse(self):

        with open(self.s_file, 'r') as fp:
            marker = False
            for line in fp.readlines():
                if '____' in line:
                    marker = True
                    continue
                if marker is False:
                    continue

                yield {'accession' : line.strip()}

    def __str__(self):
        return 'SwissWithdrawn_Parser'


class MESHChangesParser(Parser):

    def __init__(self, url):
        super(MESHChangesParser, self).__init__(url)
        self.mesh_file = url

    def parse(self):

        with open(self.mesh_file, 'r') as fp:
            for line in fp.readlines():
                if 'MH OLD =' in line:
                    mh_old = line.split('= ')[1]
                    if '#' in mh_old:
                        mh_old = mh_old.split(' #')[0]
                    elif '[' in mh_old:
                        mh_old = mh_old.split(' [')[0]
                if 'MH NEW =' in line:
                    mh_new = line.split('= ')[1]
                    if '#' in mh_new:
                        mh_new = mh_new.split(' #')[0]
                    elif '[' in mh_new:
                        mh_new = mh_new.split(' [')[0]
                    yield { 'mh_old' : mh_old.strip(), 'mh_new' : mh_new.strip() }
                    mh_old = ''
                    mh_new = ''

    def __str__(self):
        return 'MESHChanges_Parser'


def is_deprecated(child_list):
    dep = False
    deprecated = '{http://www.w3.org/2002/07/owl#}deprecated'
    for child in child_list:
        if child.tag == deprecated and child.text == 'true':
            dep = True
    return dep


# custom exception to break out of the loop when an deprecated
# term has been seen.
class DeprecatedTermException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

# excludes deprecated terms which are not included in the namespace
class DOParser(Parser):

    def __init__(self, url):
        super(DOParser, self).__init__(url)
        self.do_file = url
        self.classy = '{http://www.w3.org/2002/07/owl#}Class'
        self.deprecated = '{http://www.w3.org/2002/07/owl#}deprecated'
        self.dbxref = '{http://www.geneontology.org/formats/oboInOwl#}hasDbXref'
        self.id = '{http://www.geneontology.org/formats/oboInOwl#}id'
        self.label = '{http://www.w3.org/2000/01/rdf-schema#}label'

    def parse(self):

        with open(self.do_file, 'rb') as df:
            tree = etree.iterparse(df, tag=self.classy)
            for event, elem in tree:
                do_dict = {}
                dbxrefs = []
                name = ''
                id = ''
                try:
                    if len(elem.values()) != 0:
                        children = elem.getchildren()
                        if is_deprecated(children):
                            raise DeprecatedTermException(children)
                        else:
                            for child in children:
                                if child.tag == self.dbxref:
                                    dbxrefs.append(child.text)
                                elif child.tag == self.id:
                                    id = child.text.split(':')[1]
                                elif child.tag == self.label:
                                    name = child.text
                            do_dict['name'] = name
                            do_dict['id'] = id
                            do_dict['dbxrefs'] = dbxrefs
                            yield do_dict
                except DeprecatedTermException:
                    continue

    def __str__(self):
        return 'DO_Parser'


# includes deprecated terms (for change-log)
class DODeprecatedParser(Parser):

    def __init__(self, url):
        super(DOParser, self).__init__(url)
        self.do_file = url
        self.classy = '{http://www.w3.org/2002/07/owl#}Class'
        self.deprecated = '{http://www.w3.org/2002/07/owl#}deprecated'
        self.dbxref = '{http://www.geneontology.org/formats/oboInOwl#}hasDbXref'
        self.id = '{http://www.geneontology.org/formats/oboInOwl#}id'
        self.label = '{http://www.w3.org/2000/01/rdf-schema#}label'

    def parse(self):

        with open(self.do_file, 'rb') as df:
            tree = etree.iterparse(df, tag=self.classy)
            for event, elem in tree:
                do_dep_dict = {}
                dbxrefs = []
                name = ''
                id = ''
                dep = False
                if len(elem.values()) != 0:
                    children = elem.getchildren()
                    for child in children:
                        if child.tag == self.dbxref:
                            dbxrefs.append(child.text)
                        elif child.tag == self.id:
                            id = child.text.split(':')[1]
                        elif child.tag == self.label:
                            name = child.text
                        elif child.tag == self.deprecated:
                            dep = True
                do_dep_dict['name'] = name
                do_dep_dict['id'] = id
                do_dep_dict['dbxrefs'] = dbxrefs
                do_dep_dict['deprecated'] = dep
                yield do_dep_dict

    def __str__(self):
        return 'DODeprecated_Parser'
