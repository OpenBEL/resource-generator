# coding: utf-8

'''
 parsed.py

 Acts as a storage module for the data being parsed by each parser.
 This data can get very large, and it has been shown that this module
 alone is not sufficient to meet the memory needs of the program,
 specifically of the PubChem dataset, which is currently commented out.
 Consider replacing with DBM or another type of storage management.

'''

from datasets import *
from configuration import *

count = 0

# Data needed for the change-log
swiss_withdrawn_acc_dict = {}

swiss_withdrawn_acc_data = SwissWithdrawnData(swiss_withdrawn_acc_dict)


# entry passed to this function will be one row from
# the file being parsed by its parser.
def build_data(entry, parser, data_object):

    if parser == 'NamespaceParser':
        term_id = entry.get('ID')
        data_object._dict[term_id] = {
            'ALTIDS': entry.get('ALTIDS'),
            'LABEL': entry.get('LABEL'),
            'SYNONYMS': entry.get('SYNONYMS'),
            'DESCRIPTION': entry.get('DESCRIPTION'),
            'TYPE': entry.get('TYPE'),
            'SPECIES': entry.get('SPECIES'),
            'XREF': entry.get('XREF'),
            'OBSOLETE': entry.get('OBSOLETE'),
            'PARENTS': entry.get('PARENTS'),
            'CHILDREN': entry.get('CHILDREN')}

    if parser == 'EntrezGeneInfo_Parser':
        gene_id = entry.get('GeneID')
        type_of_gene = entry.get('type_of_gene')
        description = entry.get('description')
        tax_id = entry.get('tax_id')
        official_symbol = entry.get('Symbol_from_nomenclature_authority')
        synonyms = entry.get('Synonyms')
        alternative_names = entry.get('Other_designations')
        name = entry.get('Full_name_from_nomenclature_authority')
        symbol = entry.get('Symbol')
        description = entry.get('description')
        dbXrefs = entry.get('dbXrefs')
        data_object._dict[gene_id] = {
            'dbXrefs': dbXrefs,
            'type_of_gene': type_of_gene,
            'description': description,
            'tax_id': tax_id,
            'Symbol_from_nomenclature_authority': official_symbol,
            'Symbol': symbol,
            'Synonyms': synonyms,
            'Other_designations': alternative_names,
            'Full_name_from_nomenclature_authority': name,
            'description': description}

    elif parser == 'EntrezGeneHistory_Parser':
        gene_id = entry.get('GeneID')
        discontinued_id = entry.get('Discontinued_GeneID')
        new_id = None
        if gene_id == '-':
            status = 'withdrawn'
        else:
            status = 'retired'
            new_id = gene_id

        data_object._dict[discontinued_id] = {
            'status': status,
            'new_id': new_id
        }

    elif parser == 'Homologene_Parser':
        gene_id = entry.get('GeneID')
        tax_id = entry.get('tax_id')
        hid = entry.get('HID')
        # build dict to enable lookups by homologene group id OR gene ID
        data_object._dict['gene_ids'] = data_object._dict.get('gene_ids', {})
        data_object._dict['gene_ids'][gene_id] = {
            'homologene_group': hid,
            'tax_id': tax_id
        }
        data_object._dict['homologene_groups'] = data_object._dict.get(
            'homologene_groups', {})
        data_object._dict['homologene_groups'][
            hid] = data_object._dict.get('homologene_groups').get(hid, {})
        data_object._dict['homologene_groups'][hid][
            tax_id] = data_object._dict.get(tax_id, set())
        data_object._dict['homologene_groups'][hid][tax_id].add(gene_id)

    elif parser == 'HGNC_Parser':
        app_symb = entry.get('Approved Symbol')
        loc_type = entry.get('Locus Type')
        hgnc_id = entry.get('HGNC ID').replace('HGNC:', '')
        old_symbols = entry.get('Previous Symbols')
        old_names = entry.get('Previous Names')
        synonyms = entry.get('Synonyms')
        name_synonyms = entry.get('Name Synonyms')
        name = entry.get('Approved Name')
        mouse_ortholog = entry.get(
            'Mouse Genome Database ID (supplied by MGI)')
        rat_ortholog = entry.get('Rat Genome Database ID (supplied by RGD)')

        data_object._dict[hgnc_id] = {
            'Locus Type': loc_type,
            'Symbol': app_symb,
            'Previous Symbols': old_symbols,
            'Previous Names': old_names,
            'Synonyms': synonyms,
            'Name Synonyms': name_synonyms,
            'Approved Name': name,
            'mouse_ortholog_id': mouse_ortholog,
            'rat_ortholog_id': rat_ortholog}

    elif parser == 'MGI_Parser':
        m_symbol = entry.get('Marker Symbol')
        feature_type = entry.get('Feature Type')
        m_type = entry.get('Marker Type')
        acc_id = entry.get('MGI Accession ID').replace('MGI:', '')
        m_name = entry.get('Marker Name')
        m_syn = entry.get('Marker Synonyms (pipe-separated)')

        data_object._dict[acc_id] = {
            'Feature Type': feature_type,
            'Marker Type': m_type,
            'Symbol': m_symbol,
            'Marker Name': m_name,
            'Marker Synonyms': m_syn}

    elif parser == 'RGD_Parser':
        gene_type = entry.get('GENE_TYPE')
        name = entry.get('NAME')
        symb = entry.get('SYMBOL')
        rgd_id = entry.get('GENE_RGD_ID')
        old_symbol = entry.get('OLD_SYMBOL')
        old_name = entry.get('OLD_NAME')

        data_object._dict[rgd_id] = {
            'GENE_TYPE': gene_type,
            'NAME': name,
            'SYMBOL': symb,
            'OLD_SYMBOL': old_symbol,
            'OLD_NAME': old_name}

    elif parser == 'SwissProt_Parser':
        name = entry.get('name')
        acc = entry.get('accessions')
        primary_acc = acc[0]
        gene_type = entry.get('type')
        dbref = entry.get('dbreference')
        alt_fullnames = entry.get('alternativeFullNames')
        alt_shortnames = entry.get('alternativeShortNames')
        rec_fullname = entry.get('recommendedFullName')
        rec_shortname = entry.get('recommededShortName')
        gene_name = entry.get('geneName')
        gene_syn = entry.get('geneSynonyms')
        tax_id = entry.get('NCBI Taxonomy')

        data_object._dict[primary_acc] = {
            'name': name,
            'type': gene_type,
            'accessions': acc,
            'dbreference': dbref,
            'alternativeFullNames': alt_fullnames,
            'alternativeShortNames': alt_shortnames,
            'recommendedFullName': rec_fullname,
            'recommendedShortName': rec_shortname,
            'geneName': gene_name,
            'geneSynonym': gene_syn,
            'tax_id': tax_id}

    elif parser == 'Affy_Parser':
        probe_id = entry.get('Probe Set ID')
        entrez_gene = entry.get('Entrez Gene')
        species = entry.get('Species Scientific Name')

        data_object._dict[probe_id] = {
            'Entrez Gene': entrez_gene,
            'Species': species}

    elif parser == 'Gene2Acc_Parser':
        status = entry.get('status')
        taxid = entry.get('tax_id')
        entrez_gene = entry.get('GeneID')

        data_object._dict[entrez_gene] = {
            'status': status,
            'tax_id': taxid,
            'entrez_gene': entrez_gene}

    elif parser == 'CHEBI_Parser':
        name = entry.get('name')
        primary_id = entry.get('primary_id')
        alt_ids = entry.get('alt_ids')
        synonyms = entry.get('synonyms')

        data_object._dict[primary_id] = {
            'name': name,
            'alt_ids': alt_ids,
            'synonyms': synonyms}

    elif parser == 'GO_Parser':
        termid = entry.get('termid')
        termname = entry.get('termname')
        alt_ids = entry.get('altids')
        syn = entry.get('synonyms')
        namespace = entry.get('namespace')
        is_obsolete = entry.get('is_obsolete')
        is_complex = entry.get('complex')

        if namespace == 'cellular_component':
            gocc_dict[termid] = {
                'termname': termname,
                'alt_ids': alt_ids,
                'synonyms': syn,
                'complex': is_complex,
                'is_obsolete': is_obsolete}

        elif namespace == 'biological_process':
            gobp_dict[termid] = {
                'termname': termname,
                'alt_ids': alt_ids,
                'synonyms': syn,
                'is_obsolete': is_obsolete}

    elif parser == 'MESH_Parser':
        ui = entry.get('ui')
        mh = entry.get('mesh_header')
        mns = entry.get('mns')
        sts = entry.get('sts')
        synonyms = entry.get('synonyms')

        mesh_dict = {}
        mesh_dict['mesh_header'] = mh
        mesh_dict['synonyms'] = synonyms

        if any(branch.startswith('A11.284') for branch in mns):
            meshcl_dict[ui] = mesh_dict

        disease_branches = ('C', 'F03')
        if any(branch.startswith(disease_branches) for branch in mns):
            meshd_dict[ui] = mesh_dict

        if any(branch.startswith('G') for branch in mns):
            branches = {branch for branch in mns if branch.startswith('G')}
            excluded = ('G01', 'G15', 'G17')
            if not all(branch.startswith(excluded) for branch in branches):
                meshpp_dict[ui] = mesh_dict

        if any(branch.startswith('D') for branch in mns) or ui.startswith('C'):
            # filter by semantic type for chemicals, but allow all non-typed terms in
            # see http://semanticnetwork.nlm.nih.gov/SemGroups/SemGroups.txt
            chemicals = (
                'T116',
                'T195',
                'T123',
                'T122',
                'T118',
                'T103',
                'T120',
                'T104',
                'T200',
                'T111',
                'T196',
                'T126',
                'T131',
                'T125',
                'T129',
                'T130',
                'T197',
                'T119',
                'T124',
                'T114',
                'T109',
                'T115',
                'T121',
                'T192',
                'T110',
                'T127')

            if ui.startswith('D') or any(
                    st in chemicals for st in sts) or len(sts) == 0:
                meshc_dict[ui] = mesh_dict

        if any(branch.startswith('A') for branch in mns):
            excluded = ('A13', 'A18', 'A19', 'A20', 'A21')
            branches = {branch for branch in mns if branch.startswith('A')}
            if not all(branch.startswith(excluded) for branch in branches):
                mesha_dict[ui] = mesh_dict

    elif parser == 'SwissWithdrawn_Parser':
        term_id = entry.get('accession')

        data_object._dict[term_id] = {
            'status': 'withdrawn'
        }

    elif str(parser) == 'DO_Parser':
        name = entry.get('name')
        term_id = entry.get('id')
        dbxrefs = entry.get('dbxrefs')
        synonyms = entry.get('synonyms')
        data_object._dict[term_id] = {
            'name': name,
            'dbxrefs': dbxrefs,
            'synonyms': synonyms,
            'alt_ids': entry.get('alt_ids'),
            'is_obsolete': entry.get('is_obsolete')}

    elif parser == 'Owl_Parser':
        term_id = entry.get('id')
        prefix = data_object._prefix.upper()
        term_type = entry.get('term_type')
        if prefix == 'EFO' and term_type is None:
            pass
        elif term_id.startswith(prefix):
            term_id = term_id.lstrip(prefix + '_')
            term_id = term_id.lstrip(prefix + 'ID_')
            data_object._dict[term_id] = {
                'name': entry.get('name'),
                'dbxrefs': entry.get('dbxrefs'),
                'synonyms': entry.get('synonyms'),
                'alt_ids': entry.get('alt_ids'),
                'is_obsolete': entry.get('is_obsolete')}
            if term_type:
                data_object._dict[term_id]['term_type'] = term_type

    elif parser == 'RGDOrthologParser':
        term_id = entry.get('RAT_GENE_RGD_ID')
        data_object._dict[term_id] = {
            'symbol': entry.get('RAT_GENE_SYMBOL'),
            'human_ortholog_id': entry.get('HUMAN_ORTHOLOG_HGNC_ID'),
            'mouse_ortholog_id': entry.get('MOUSE_ORTHOLOG_MGI')}

    elif parser == 'RGD_Obsolete_Parser':
        term_id = entry.get('OLD_GENE_RGD_ID')
        data_object._dict[term_id] = {
            'status': entry.get('OLD_GENE_STATUS').lower(),
            'new_id': entry.get('NEW_GENE_RGD_ID'),
            'type': entry.get('OLD_GENE_TYPE')	}

    elif parser == 'NCBI_Taxonomy_Parser':
        term_id = entry.get('term_id')
        data_object._dict[term_id] = {
            'name': entry.get('name'),
            'synonyms': entry.get('synonyms')}

