from string import Template

# templates for annotation blocks
author = Template("""[Author]
NameString=OpenBEL
CopyrightString=Copyright (c) ${year}, OpenBEL Project. This work is licensed under a Creative Commons Attribution 3.0 Unported License.
ContactInfoString=info@openbel.org\n
""")

processing = """[Processing]
CaseSensitiveFlag=yes
DelimiterString=|
CacheableFlag=yes\n\n"""


def annotation_definition(name, version, crdate):
    """ Generate [AnnotationDefinition] block, given annotation name,
    .belanno version, and creation date. """
    annodef_template = Template("""[AnnotationDefinition]
Keyword=${kw}
TypeString=list
DescriptionString=${ad}
UsageString=${us}
VersionString=${version}
CreatedDateTime=$date\n
""")
    hd = header_dict[name]
    annodef = annodef_template.substitute(
        kw=hd['kw'],
        ad=hd['ad'],
        us=hd['us'],
        version=version,
        date=crdate,)
    return annodef


def citation_info(id, source_version, source_date=None):
    """ Generate information for [Citation] block, given
    source id, version, and date (optional). """
    data = source_data[id]
    citation = Template("""NameString=${cname}
DescriptionString=${sd}
PublishedVersionString=${sv}\n""").substitute(
        cname=data['cname'],
        sd=data['sd'],
        sv=source_version)
    if source_date:
        citation += Template("""PublishedDate=${sdate}\n""").substitute(
            sdate=source_date)
    citation += Template("""ReferenceURL=${su}\n""").substitute(
        su=data['su'])
    return citation

# header_dict contains information for the [AnnotationDefinition] block
header_dict = {}
header_dict['anatomy'] = {
    'kw': 'Anatomy',
    'ad': 'Anatomy terms from Uberon',
    'us': 'Use to annotate a statement demonstrated in a specific anatomical entity.'
}

header_dict['disease'] = {
    'kw': 'Disease',
    'ad': 'Disease names from the Disease Ontology (DO)',
    'us': 'Use to annotate a statement demonstrated in a biological system where the disease state was present or where the facts represented in the statement were shown to be relevant to the disease.'
}

header_dict['cell'] = {
    'kw': 'Cell',
    'ad': 'Cell types from the Cell Ontology (CL)',
    'us': 'Use this annotation to indicate the cell type context of a statement. Use to annotate a statement that was demonstrated in a given cell type or where the facts represented in the statement were shown to be relevant to the cell type.'
}

header_dict['cell-line'] = {
    'kw': 'CellLine',
    'ad': 'Cell Line Ontology (CLO) and Experimental Factor Ontology (EFO) Cell Lines',
    'us': 'Use this annotation to indicate the cell line context of a statement. Use to annotate a statement that was demonstrated in a given cell line.'
}

header_dict['cell-structure'] = {
    'kw': 'CellStructure',
    'ad': 'The Cell Structures [A11.284] sub-branch of the Medical Subject Headings (MeSH) Anatomy [A] branch.',
    'us': 'Use to annotate a statement demonstrated in the context of a specific cellular location.'}

header_dict['mesh-diseases'] = {
    'kw': 'MeSHDisease',
    'ad': 'Disease terms from the [C] and [F03] branches of Medical Subject Headings (MeSH).',
    'us': 'Use to annotate a statement demonstrated in the context of a specific disease.',
}

header_dict['mesh-anatomy'] = {
    'kw': 'MeSHAnatomy',
    'ad': 'Anatomy terms from Medical Subject Headings (MeSH), limited to headings from the following subbranches:[A01], [A02], [A03], [A04], [A05], [A06], [A07], [A08], [A09], [A10], [A12], [A14], [A15], [A16], and [A17].',
    'us': 'Use to annotate a statement demonstrated in a specific anatomical entity.',
    'cname': '',
    'sd': '',
    'su': ''
}

# source_data contains information for the [Citation] block
source_data = {}
source_data['UBERON'] = {
    'cname': 'Uberon',
    'sd': 'Uberon is an integrated cross-species anatomy ontology representing a variety of entities classified according to traditional anatomical criteria such as structure, function and developmental lineage.',
    'su': 'http://uberon.org/'
}

source_data['DOID'] = {
    'cname': 'Disease Ontology (DO)',
    'sd': 'The Disease Ontology has been developed as a standardized ontology for human disease with the purpose of providing the biomedical community with consistent, reusable and sustainable descriptions of human disease terms, phenotype characteristics and related medical vocabulary disease concepts through collaborative efforts of researchers at Northwestern University, Center for Genetic Medicine and the University of Maryland School of Medicine, Institute for Genome Sciences. The Disease Ontology semantically integrates disease and medical vocabularies through extensive cross mapping of DO terms to MeSH, ICD, NCI\'s thesaurus, SNOMED and OMIM.',
    'su': 'http://disease-ontology.org/'
}

source_data['CL'] = {
    'cname': 'Cell Ontology (CL)',
    'sd': 'The Cell Ontology (CL) is a candidate OBO Foundry ontology for the representation of cell types. First described in 2005, the CL integrates cell types from the prokaryotic, fungal, and eukaryotic organisms.',
    'su': 'http://cellontology.org/'
}

source_data['CLO'] = {
    'cname': 'Cell Line Ontology (CLO)',
    'sd': 'The Cell Line Ontology (CLO) is a community-driven ontology that is developed to standardize and integrate cell line information and support computer-assisted reasoning.',
    'su': 'http://www.clo-ontology.org/'}

source_data['EFO'] = {
    'cname': 'Experimental Factor Ontology (EFO)',
    'sd': 'The Experimental Factor Ontology (EFO) provides a systematic description of many experimental variables available in EBI databases.',
    'su': 'http://www.ebi.ac.uk/efo/'}

source_data['MESH'] = {
    'cname': 'MeSH',
    'sd': 'MeSH (Medical Subject Headings) is a controlled vocabulary thesaurus created, maintained, and provided by the U.S. National Library of Medicine.',
    'su': 'http://www.nlm.nih.gov/mesh/meshhome.html'}

