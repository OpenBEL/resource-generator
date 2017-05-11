# coding: utf-8

'''
 configuration.py

 Data objects are initialized here; these objects are expected to be named
 using the data object '_prefix' string attribute + '_data'.

 baseline_data - ordered dictionary mapping each source/input file
 to its parser and resulting data objects. The local file name is
 mapped to a data tuple containing:
    [0] the source file url
    [1] the parser
    [2] the data object, or a list of data objects
        in the case that multiple objects are generated from the same file
 This dictionary is consumed by gp_baseline.

 data_file_info - dictionary mapping source file names
 to .belns/.beleq output files. This is used by common.get_citation_info
 only in the case where a different source file than the one used to generate
 the corresponding data object should be cited, e.g, the source of the HGNC/RGD/MGI
 equivalence mappings used for the HGNC/MGI/RGD .beleq files is EntrezGene.

 affy_array_names - list of array names to include in the AFFX (Affymetrix
 probe set) namespace. This list is consumed by parsers.AffyParser.


'''

from collections import OrderedDict
from common import get_latest_GO_filename
from common import get_latest_MeSH_filename
from datasets import *
import parsers
import os

file_url = 'file://{0}/datasets/'.format(os.getcwd())

baseline_data = OrderedDict()

egid_data = EntrezInfoData()
baseline_data['entrez_info.gz'] = (
    'ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/All_Mammalia.gene_info.gz',
    parsers.EntrezGeneInfoParser,
    egid_data)

entrez_history_data = EntrezHistoryData()
baseline_data['entrez_history.gz'] = (
    'ftp://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz',
    parsers.EntrezGeneHistoryParser, entrez_history_data
)

hgnc_data = HGNCData()
baseline_data['hgnc.tsv'] = (
    'ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc_complete_set.txt.gz',
    parsers.HGNCParser, hgnc_data
)

mgi_data = MGIData()
baseline_data['mgi.rpt'] = (
    'http://www.informatics.jax.org/downloads/reports/MRK_List2.rpt',
    parsers.MGIParser, mgi_data
)

# need to create classes MGIHistoryData, MGIHistoryParser
# mgi_history_data = MGIHistoryData()
# baseline_data[mgi_history_index] = (
# 'ftp://ftp.informatics.jax.org/pub/reports/archive/nomen/',
# parsers.MGIHistoryParser, mgi_history_data
# )

rgd_data = RGDData()
baseline_data['rgd.txt'] = (
    'ftp://ftp.rgd.mcw.edu/pub/data_release/GENES_RAT.txt',
    parsers.RGDParser, rgd_data
)

rgd_history_data = RGDObsoleteData()
baseline_data['rgd_obsolete.txt'] = (
    'ftp://ftp.rgd.mcw.edu/pub/data_release/GENES_OBSOLETE_IDS.txt',
    parsers.RGDObsoleteParser,
    rgd_history_data)

rgd_ortho_data = OrthologyData(prefix='rgd')
baseline_data['RGD_ORTHOLOGS.txt'] = (
    'ftp://ftp.rgd.mcw.edu/pub/data_release/RGD_ORTHOLOGS.txt',
    parsers.RGDOrthologParser,
    rgd_ortho_data)

sp_data = SwissProtData()
baseline_data['swiss.xml.gz'] = (
    'ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz',
    parsers.SwissProtParser,
    sp_data)

sp_history_data = HistoryDataSet(prefix='sp')
baseline_data['delac_sp.txt'] = (
    'ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/docs/delac_sp.txt',
    parsers.SwissWithdrawnParser,
    sp_history_data)

gene2acc_data = Gene2AccData()
baseline_data['gene2acc.gz'] = (
    'ftp://ftp.ncbi.nih.gov/gene/DATA/gene2accession.gz',
    parsers.Gene2AccParser, gene2acc_data
)


affx_data = AffyData()
baseline_data['affy.xml'] = (
    'http://www.affymetrix.com/analysis/downloads/netaffxapi/GetFileList.jsp?licence=OPENBEL2013&user=jhourani@selventa.com&password=OPENBEL2013',
    parsers.AffyParser,
    affx_data)
chebi_data = CHEBIData()
baseline_data['chebi.owl'] = (
    'ftp://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.owl',
    parsers.CHEBIParser, chebi_data
)

schem_data = StandardCustomData(
    name='selventa-legacy-chemicals',
    prefix='schem',
    domain=['chemical'])
baseline_data['selventa-legacy-chemical-names.txt'] = (
    file_url +
    'selventa-legacy-chemical-names.txt',
    parsers.NamespaceParser,
    schem_data)

sdis_data = StandardCustomData(
    name='selventa-legacy-diseases',
    prefix='sdis',
    domain=['chemical'])
baseline_data['selventa-legacy-diseases.txt'] = (
    file_url +
    'selventa-legacy-diseases.txt',
    parsers.NamespaceParser,
    sdis_data)

scomp_data = StandardCustomData(
    name='selventa-named-complexes',
    prefix='scomp',
    domain=['complex'])
baseline_data['selventa-named-complexes.txt'] = (
    file_url +
    'selventa-named-complexes.txt',
    parsers.NamespaceParser,
    scomp_data)

sfam_data = StandardCustomData(
    name='selventa-protein-families',
    prefix='sfam',
    domain=['other'])
baseline_data['selventa-protein-families.txt'] = (
    file_url +
    'selventa-protein-families.txt',
    parsers.NamespaceParser,
    sfam_data)

tax_data = NCBITaxonomyData(
    name='ncbi-taxonomy',
    prefix='taxon',
    domain=['species'],
    scheme_type=['anno'])
baseline_data['taxnames_9606_10090_10116.dmp'] = (
    file_url +
    'taxnames_9606_10090_10116.dmp',
    parsers.NCBITaxonomyParser,
    tax_data)

# - get the latest GO archive file name and URL
#go_file = get_latest_GO_filename('http://archive.geneontology.org/latest-full')
go_file = 'http://archive.geneontology.org/latest-full/go_monthly-termdb.obo-xml.gz'
gobp_dict, gocc_dict = {}, {}
gobp_data = GOData(
    gobp_dict,
    name='go-biological-process',
    prefix='gobp',
    domain=['process'])
gocc_data = GOData(
    gocc_dict,
    name='go-cellular-component',
    prefix='gocc',
    domain=[
        'complex',
         'location'])
baseline_data['go.xml.gz'] = (
    go_file, parsers.GOParser, [
        gobp_data, gocc_data])

mesh_file = get_latest_MeSH_filename(
    'ftp://nlmpubs.nlm.nih.gov/online/mesh/.asciimesh/', 'd', '.bin')
meshcl_dict, meshd_dict, meshpp_dict, meshc_dict, mesha_dict = {}, {}, {}, {}, {}
meshcl_data = MESHData(
    meshcl_dict,
    name='mesh-cellular-structures',
    prefix='meshcs',
    domain=['location'],
    scheme_type=[
        'ns',
         'anno'])
meshd_data = MESHData(
    meshd_dict,
    name='mesh-diseases',
    prefix='meshd',
    domain=['disease'],
    scheme_type=[
        'ns',
         'anno'])
meshpp_data = MESHData(
    meshpp_dict,
    name='mesh-processes',
    prefix='meshpp',
    domain=['process'])
meshc_data = MESHData(
    meshc_dict,
    name='mesh-chemicals',
    prefix='meshc',
    domain=['chemical'])
mesha_data = MESHData(
    mesha_dict,
    name='mesh-anatomy',
    prefix='mesha',
    domain=['anatomy'],
    scheme_type=['anno'])
baseline_data['mesh.bin'] = (
    mesh_file, parsers.MESHParser, [
        meshcl_data, meshd_data, meshpp_data, meshc_data, mesha_data])
mesh_sup_file = get_latest_MeSH_filename(
    'ftp://nlmpubs.nlm.nih.gov/online/mesh/.asciimesh/', 'c2', '.bin')
baseline_data['meshc.bin'] = (
    mesh_sup_file,
    parsers.MESHParser, [meshc_data]
)

do_data = OWLData(
    name='disease-ontology',
    prefix='do',
    domain=['disease'],
    scheme_type=[
        'ns',
         'anno'])
baseline_data['doid.owl'] = (
    'http://purl.obolibrary.org/obo/doid.owl',
    parsers.OwlParser, do_data
)

efo_data = OWLData(
    prefix='efo',
    name='experimental-factor-ontology',
    domain=['cell-line'],
    scheme_type=['anno'])
baseline_data['efo.owl'] = (
    'http://www.ebi.ac.uk/efo/efo.owl',
    parsers.OwlParser,
    efo_data)

cl_data = OWLData(
    name='cell-ontology',
    prefix='cl',
    domain=['cell'],
    scheme_type=['anno'])
baseline_data['cl.owl'] = (
    'http://purl.obolibrary.org/obo/cl.owl',
    parsers.OwlParser,
    cl_data)

clo_data = OWLData(
    prefix='clo',
    name='cell-line-ontology',
    domain=['cell-line'],
    scheme_type=['anno'])
baseline_data['clo.owl'] = (
    'http://purl.obolibrary.org/obo/clo.owl',
    parsers.OwlParser,
    clo_data)

uberon_data = OWLData(
    name='uberon',
    prefix='uberon',
    domain=['anatomy'],
    scheme_type=['anno'])
baseline_data['uberon.owl'] = (
    'http://purl.obolibrary.org/obo/uberon.owl',
    parsers.OwlParser,
    uberon_data)
###
homologene_ortho_data = HomologeneData()
baseline_data['homologene.data'] = (
    'ftp://ftp.ncbi.nlm.nih.gov/pub/HomoloGene/current/homologene.data',
    parsers.HomologeneParser,
    homologene_ortho_data
)

affy_array_names = [
    'HG-U133A',
    'HG-U133B',
    'HG-U133_Plus_2',
    'HG_U95Av2',
    'MG_U74A',
    'MG_U74B',
    'MG_U74C',
    'MOE430A',
    'MOE430B',
    'Mouse430A_2',
    'Mouse430_2',
    'RAE230A',
    'RAE230B',
    'Rat230_2',
    'HT_HG-U133_Plus_PM',
    'HT_MG-430A',
    'HT_Rat230_PM',
    'MG_U74Av2',
    'MG_U74Bv2',
    'MG_U74Cv2',
    'HG-U219']

# location of source file to be associated with each .belns/.beleq file
#  - ONLY if different from the source file associated with the namespace data;
#  used for common.get_citation_info
data_file_info = {
    'hgnc-human-genes.beleq': 'entrez_info.gz',
    'mesh-diseases.beleq': 'doid.owl',
    'mesh-diseases-ids.beleq': 'doid.owl',
    'mgi-mouse-genes.beleq': 'entrez_info.gz',
    'rgd-rat-genes.beleq': 'entrez_info.gz',
}

