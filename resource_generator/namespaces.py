# coding: utf-8
#
# namespaces.py

import ipdb

# namespace dictionaries
entrez_ns = set()
hgnc_ns = set()
mgi_ns  = set()
rgd_ns = set()
sp_ns = set()
sp_acc_ns = set()
affy_ns = set()
chebi_name_ns = set()
chebi_id_ns = set()
pub_ns = set()

# miscRNA should not be used here, as it will be handled in a special case.
# For completion sake it is included.
entrez_encoding = {'protein-coding' : 'GRP', 'miscRNA' : 'GR', 'ncRNA' : 'GR',
                   'snoRNA' : 'GR', 'snRNA' : 'GR', 'tRNA' : 'GR',
                   'scRNA' : 'GR', 'other' : 'G', 'pseudo' : 'GR',
                   'unknown' : 'GRP', 'rRNA' : 'GR'}

hgnc_encoding = {'gene with protein product' : 'GRP', 'RNA, cluster' : 'GR',
                 'RNA, long non-coding' : 'GR', 'RNA, micro' : 'GRM',
                 'RNA, ribosomal' : 'GR', 'RNA, small cytoplasmic' : 'GR',
                 'RNA, small misc' : 'GR', 'RNA, small nuclear' : 'GR',
                 'RNA, small nucleolar' : 'GR', 'RNA, transfer' : 'GR',
                 'phenotype only' : 'G', 'RNA, pseudogene' : 'GR',
                 'T cell receptor pseudogene' : 'GR',
                 'immunoglobulin pseudogene' : 'GR', 'pseudogene' : 'GR',
                 'T cell receptor gene' : 'GRP',
                 'complex locus constituent' : 'GRP',
                 'endogenous retrovirus' : 'G', 'fragile site' : 'G',
                 'immunoglobulin gene' : 'GRP', 'protocadherin' : 'GRP',
                 'readthrough' : 'GR', 'region' : 'G',
                 'transposable element' : 'G', 'unknown' : 'GRP',
                 'virus integration site' : 'G', 'RNA, micro' : 'GRM',
                 'RNA, misc' : 'GR', 'RNA, Y' : 'GR', 'RNA, vault' : 'GR',
                 }

mgi_encoding = {'gene' : 'GRP', 'protein coding gene' : 'GRP',
                'non-coding RNA gene' : 'GR', 'rRNA gene' : 'GR',
                'tRNA gene' : 'GR', 'snRNA gene' : 'GR', 'snoRNA gene' : 'GR',
                'miRNA gene' : 'GRM', 'scRNA gene' : 'GR',
                'lincRNA gene' : 'GR', 'RNase P RNA gene' : 'GR',
                'RNase MRP RNA gene' : 'GR', 'telomerase RNA gene' : 'GR',
                'unclassified non-coding RNA gene' : 'GR',
                'heritable phenotypic marker' : 'G', 'gene segment' : 'G',
                'unclassified gene' : 'GRP', 'other feature types' : 'G',
                'pseudogene' : 'GR', 'transgene' : 'G',
                'other genome feature' : 'G', 'pseudogenic region' : 'GR',
                'polymorphic pseudogene' : 'GRP',
                'pseudogenic gene segment' : 'GR', 'SRP RNA gene' : 'GR'}

rgd_encoding = {'gene' : 'GRP', 'miscrna' : 'GR', 'predicted-high' : 'GRP',
                'predicted-low' : 'GRP', 'predicted-moderate' : 'GRP',
                'protein-coding' : 'GRP', 'pseudo' : 'GR', 'snrna' : 'GR',
                'trna' : 'GR', 'rrna' : 'GR'}

hgnc_map = {}
mgi_map = {}
rgd_map = {}
# takes a dataset 'Object' and build namespace
def make_namespace(d):

    # build and write out the namespace values
    delim = '|'
    if str(d) == 'entrez_info':
        with open('entrez-info_namespace.belns', 'w') as fp:
            # tuple of (gene_id, gene_type, description)
            for vals in d.get_ns_values():
                gene_id, gene_type, description = vals
                if gene_type == 'miscRNA':
                    if 'microRNA' in description:
                        fp.write(delim.join((gene_id, 'GRM'))+'\n')
                        entrez_ns.add(gene_id)
                    else:
                        fp.write(delim.join((gene_id, 'GR'))+'\n')
                        entrez_ns.add(gene_id)
                else:
                    fp.write(delim.join((gene_id, entrez_encoding[gene_type]))+'\n')
                    entrez_ns.add(gene_id)

    elif str(d) == 'hgnc':
        with open('hgnc-namespace.belns', 'w') as fp:
            for vals in d.get_ns_values():
                approved_symb, locus_type, hgnc_id = vals
                # withdrawn genes NOT included in this namespace
                if locus_type is not 'withdrawn' and 'withdrawn' not in approved_symb:
                    fp.write(delim.join((approved_symb, hgnc_encoding[locus_type]))+'\n')
                    hgnc_ns.add(approved_symb)
                hgnc_map[hgnc_id] = approved_symb

    elif str(d) == 'mgi':
        with open('mgi-namespace.belns', 'w') as fp:
            for vals in d.get_ns_values():
                marker_symbol, feature_type, acc_id, marker_type = vals
                if marker_type == 'Gene' or marker_type == 'Pseudogene':
                    fp.write(delim.join((marker_symbol, mgi_encoding[feature_type]))+'\n')
                    mgi_ns.add(marker_symbol)
                mgi_map[acc_id] = marker_symbol

    # withdrawn genes are NOT included in this namespace
    elif str(d) == 'rgd':
        with open('rgd-namespace.belns', 'w') as fp:
            for vals in d.get_ns_values():
                symbol, gene_type, name, rgd_id = vals
                if gene_type == 'miscrna' and 'microRNA' in name:
                    fp.write(delim.join((symbol, 'GRM'))+'\n')
                    rgd_ns.add(symbol)
                elif gene_type == 'miscrna' and 'microRNA' not in name:
                    fp.write(delim.join((symbol, 'GR'))+'\n')
                    rgd_ns.add(symbol)
                else:
                    if gene_type is not '':
                        fp.write(delim.join((symbol, rgd_encoding[gene_type]))+'\n')
                        rgd_ns.add(symbol)
                rgd_map[rgd_id] = symbol

    elif str(d) == 'swiss':
        with open('swiss-namespace.belns', 'w') as fp, open('swiss-acc-namespace.belns', 'w') as f:
            for vals in d.get_ns_values():
                gene_name, accessions = vals
                fp.write(delim.join((gene_name, 'GRP'))+'\n')
                sp_ns.add(gene_name)
                for acc in accessions:
                    f.write(delim.join((acc, 'GRP'))+'\n')
                    sp_acc_ns.add(acc)

    # are there duplicates being taken in here??
    elif str(d) == 'affy':
        with open('affy-namespace.belns', 'w') as fp:
            for vals in d.get_ns_values():
                probe_set_ids = vals
                for pid in probe_set_ids:
                    fp.write(delim.join((pid, 'R'))+'\n')
#                    if pid not in affy_ns_dict:
#                        affy_ns_dict[pid] = 'R'

    elif str(d) == 'chebi':
        with open('chebi-namespace.belns', 'w') as fp, open('chebi-id-namespace.belns', 'w') as f:
            for vals in d.get_ns_values():
                name, primary_id, altIds = vals
                fp.write(delim.join((name, 'A'))+'\n')
                chebi_name_ns.add(name)
                f.write(delim.join((primary_id, 'A'))+'\n')
                chebi_id_ns.add(name)
                if altIds:
                    for i in altIds:
                        # this check should not be needed...(more pythonic way?)
                        if i not in chebi_id_ns:
                            f.write(delim.join((i, 'A'))+'\n')
                        chebi_id_ns.add(i)

    elif str(d) == 'pubchem_namespace':
        with open('pubchem-namespace.belns', 'w') as fp:
            for vals in d.get_ns_values():
                pid = vals
                fp.write(delim.join((pid, 'A'))+'\n')
                pub_ns.add(pid)
