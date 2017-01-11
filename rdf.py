#!/usr/bin/env python3
# coding: utf-8

import argparse
import os
import pickle
import datasets
from collections import defaultdict
from rdflib import URIRef, Literal, Namespace, Graph
from rdflib.namespace import RDF, SKOS, DCTERMS, XSD
from urllib import parse

namespace = Namespace("http://www.openbel.org/bel/namespace/")
BELV = Namespace("http://www.openbel.org/vocabulary/")

annotation_concept_types = {
    'AnnotationConcept',
    'AnatomyAnnotationConcept',
    'CellLineAnnotationConcept',
    'CellAnnotationConcept',
    'DiseaseAnnotationConcept',
    'LocationAnnotationConcept',
    'SpeciesAnnotationConcept'}


def literal(obj):
    if isinstance(obj, str):
        return Literal(str(obj), datatype=XSD.string)
    else:
        return Literal(obj)

# loads parsed data from pickle file (after running phase 2 of gp_baseline.py)


def make_rdf(d, g, prefix_dict=None):
    ''' Given namepsace data object d and graph g,
    adds to  namespace rdf graph. '''
    # make namespace for data set (using data object attributes '_name' and
    # '_prefix')
    n = Namespace("http://www.openbel.org/bel/namespace/" + d._name + '/')

    print('building RDF graph for {0} ...'.format(n))
    # bind namespace prefixes
    g.bind("skos", SKOS)
    g.bind("dcterms", DCTERMS)
    g.bind("belv", BELV)
    g.bind(d._prefix, n)

    if 'ns' in d.scheme_type:
        g.add((namespace[d._name], RDF.type, BELV.NamespaceConceptScheme))
    if 'anno' in d.scheme_type:
        g.add((namespace[d._name], RDF.type, BELV.AnnotationConceptScheme))

    name = d._name.replace('-', ' ').title()
    domain = d._domain
    g.add((namespace[d._name], SKOS.prefLabel, literal(name)))
    g.add((namespace[d._name], BELV.prefix, literal(d._prefix.upper())))
    # TODO consider updating domain to uri instead of literal
    # domain is a list - a namespace can be associated with multiple domains
    for dom in domain:
        g.add((namespace[d._name], BELV.domain, literal(dom)))

    for term_id in d.get_values():
        term_clean = parse.quote(term_id)
        term_uri = URIRef(n[term_clean])
        # add primary identifier
        g.add((term_uri, DCTERMS.identifier, literal(term_id)))
        # add secondary/alternative identifiers (alt_ids)
        alt_ids = d.get_alt_ids(term_id)
        if alt_ids:
            for alt_id in alt_ids:
                g.add((term_uri, DCTERMS.identifier, literal(alt_id)))
                # adding history info (gives each alt_id a URI and links
                # to primary id URI)
                alt_id_clean = parse.quote(alt_id)
                alt_uri = URIRef(n[alt_id_clean])
                g.add((term_uri, DCTERMS.replaces, alt_uri))
                g.add((alt_uri, BELV.status, literal('secondary')))
                g.add((alt_uri, SKOS.inScheme, namespace[d._name]))

        # add official name (as title)
        name = d.get_name(term_id)
        if name:
            g.add((term_uri, DCTERMS.title, literal(name)))

        # link to to Concept Scheme
        g.add((term_uri, SKOS.inScheme, namespace[d._name]))
        pref_label = d.get_label(term_id)
        if pref_label:
            g.add((term_uri, SKOS.prefLabel, literal(pref_label)))

        # add species (tax id as literal)
        species = d.get_species(term_id)
        if species:
            g.add((term_uri, BELV.fromSpecies, literal(species)))

        # use encoding information to determine concept types
        if 'ns' in d.scheme_type:
            encoding = d.get_encoding(term_id)
            if 'G' in encoding:
                g.add((term_uri, RDF.type, BELV.GeneConcept))
            if 'R' in encoding:
                g.add((term_uri, RDF.type, BELV.RNAConcept))
            if 'M' in encoding:
                g.add((term_uri, RDF.type, BELV.MicroRNAConcept))
            if 'P' in encoding:
                g.add((term_uri, RDF.type, BELV.ProteinConcept))
            if 'A' in encoding:
                g.add((term_uri, RDF.type, BELV.AbundanceConcept))
            if 'B' in encoding:
                g.add((term_uri, RDF.type, BELV.BiologicalProcessConcept))
            if 'C' in encoding:
                g.add((term_uri, RDF.type, BELV.ComplexConcept))
            if 'O' in encoding:
                g.add((term_uri, RDF.type, BELV.PathologyConcept))
        # add conceptType (for now, just annotations)
        concept_types = d.get_concept_type(term_id)

        if 'anno' in d.scheme_type:
            concept_types = {c + 'AnnotationConcept' for c in concept_types}
            for c in concept_types:
                if c in annotation_concept_types:
                    g.add((term_uri, RDF.type, BELV[c]))
                else:
                    print(
                        'WARNING - {0} not a BELV AnnotationConcept'.format(c))

        # get synonyms (alternative symbols and names)
        alt_symbols = d.get_alt_symbols(term_id)
        if alt_symbols:
            for symbol in alt_symbols:
                g.add((term_uri, SKOS.altLabel, literal(symbol)))
        alt_names = d.get_alt_names(term_id)
        if alt_names:
            for name in alt_names:
                g.add((n[term_id], SKOS.altLabel, literal(name)))

        # get equivalences to other namespaces (must be in data set)
        if prefix_dict is None:
            print('Building required prefix dictionary ...')
            prefix_dict = build_prefix_dict()
        xrefs = d.get_xrefs(term_id)
        if xrefs:
            xref_dict = defaultdict(set)
            for x in xrefs:
                # xrefs are expected in format PREFIX:value
                if len(x.split(':')) == 2:
                    prefix, value = x.split(':')
                    xref_dict[prefix].add(value)
            for prefix, values in xref_dict.items():
                # NOTE - only xrefs with prefixes that match namespaces in the data set will be used
                # if term has multiple xrefs in same namespace, use closeMatch
                # instead of exactMatch!
                if prefix.lower() in prefix_dict:
                    xrefns = Namespace(
                        "http://www.openbel.org/bel/namespace/" +
                        prefix_dict[
                            prefix.lower()] +
                        '/')
                    g.bind(prefix.lower(), xrefns)
                    if len(values) == 1:
                        xref_uri = URIRef(xrefns[values.pop()])
                        g.add((term_uri, SKOS.exactMatch, xref_uri))
                    elif len(values) > 1:
                        for value in values:
                            xref_uri = URIRef(xrefns[value])
                            g.add((term_uri, SKOS.closeMatch, xref_uri))


def get_close_matches(concept_type, g):
    ''' Given namespace concept type and graph, g,
    searches graph for items of the given concept type and
    identifies close matches based on case-insensitive string matching
    of synonyms. Adds skos:closeMatch edges to graph. '''

    print('gathering string matches for {0}s ...'.format(concept_type))
    concept_dict = defaultdict(dict)
    for c in g.subjects(RDF.type, BELV[concept_type]):
        labels = list(g.objects(c, SKOS["altLabel"]))
        labels.extend(g.objects(c, SKOS["prefLabel"]))
        labels = [l.casefold() for l in labels]
        scheme = set(g.objects(c, SKOS["inScheme"])).pop()
        concept_dict[scheme][c] = set(labels)

    count = 0
    for scheme, concept_map in concept_dict.items():
        domains = set(g.objects(scheme, BELV["domain"]))
        for concept, labels in concept_map.items():
            # iterate concept_dict a 2nd time, skipping items in same scheme
            for scheme2 in concept_dict.keys():
                domains2 = set(g.objects(scheme2, BELV["domain"]))
                if scheme2 == scheme:
                    continue
                if domains.isdisjoint(
                        domains2):  # skip if no domain matches for scheme and scheme2
                    continue
                else:
                    for concept2, labels2 in concept_dict[scheme2].items():
                        if len(labels.intersection(labels2)) > 0:
                            g.add((concept, SKOS.closeMatch, concept2))
                            count += 1
    print('\tadded {0} closeMatches for {1}s'.format(count, concept_type))


def get_ortho_matches(d, g, prefix_dict=None):
    ''' Given an OrthologyData object and graph, adds orthologousMatch
    relationships to graph. Prefixes for both the data object and
    orthologs must be in the prefix dictionary. If prefix_dict is
    not provided, it will be generated from the namespace data objects
    in the working directory. '''

    print('gathering orthology information from {0} ...'.format(str(d)))
    if not isinstance(d, datasets.OrthologyData):
        print('ERROR - {0} is not an OrthologyData object!'.format(str(d)))
        return None
    else:
        if prefix_dict is None:
            print('Building required prefix dictionary ...')
            prefix_dict = build_prefix_dict()
        n = Namespace(
            "http://www.openbel.org/bel/namespace/" +
            prefix_dict[
                d._prefix] +
            '/')
        # for term_id in d._dict.keys():
        for term_id in d.get_values():
            term_clean = parse.quote(term_id)
            term_uri = URIRef(n[term_clean])
            ortho = d.get_orthologs(term_id)
            if ortho is not None:
                for o in ortho:
                    if len(o.split(':')) == 2:
                        prefix, value = o.split(':')
                        if prefix.lower() in prefix_dict:
                            ortho_ns = Namespace(
                                "http://www.openbel.org/bel/namespace/" +
                                prefix_dict[
                                    prefix.lower()] +
                                '/')
                            g.bind(prefix.lower(), ortho_ns)
                            ortho_uri = URIRef(ortho_ns[value])
                            g.add((term_uri, BELV.orthologousMatch, ortho_uri))


def get_history_data(d, g, prefix_dict=None):
    ''' Given a HistoryData object and graph, add status and replaces replationships
     to graph.'''

    print('gathering history information from {0} ...'.format(str(d._prefix)))
    if not isinstance(d, datasets.HistoryDataSet):
        print('ERROR - {0} is not a HistoryDataSet object!'.format(str(d)))
        return None

    else:
        if prefix_dict is None:
            print('Building required prefix dictionsry ...')
            prefix_dict = build_prefix_dict()
        n = Namespace(
            "http://www.openbel.org/bel/namespace/" +
            prefix_dict[
                d._prefix] +
            '/')
        scheme = namespace[prefix_dict[d._prefix]]
        obsolete = d.get_obsolete_ids()
        for term_id, new_value in obsolete.items():
            term_clean = parse.quote(term_id)
            term_uri = URIRef(n[term_clean])
            g.add((term_uri, SKOS.inScheme, scheme))
            if new_value == 'withdrawn':
                g.add((term_uri, BELV.status, literal('withdrawn')))
            elif new_value is not None:
                new_value_clean = parse.quote(new_value)
                new_uri = URIRef(n[new_value_clean])
                g.add((term_uri, BELV.status, literal('retired')))
                g.add((new_uri, DCTERMS.replaces, term_uri))
            else:
                print('Check values for {0}: {1}'.format(str(d), term_id))


def build_prefix_dict():
    ''' Build dictionary of namespace prefixes to names (from data set objects). '''
    print('gathering namespace information ...')
    prefix_dict = {}
    for files in os.listdir("."):
        if files.endswith("parsed_data.pickle"):
            with open(files, 'rb') as f:
                d = pickle.load(f)
            if isinstance(d, datasets.NamespaceDataSet):
                prefix_dict[d._prefix] = d._name
                print('\t{0} - {1}'.format(d._prefix, d._name))
    return prefix_dict

if __name__ == '__main__':
    # command line arguments - directory for pickled data objects
    parser = argparse.ArgumentParser(description="""Generate namespace and equivalence files
	for gene/protein datasets.""")

    parser.add_argument(
        "-n",
        required=True,
        metavar="DIRECTORY",
        help="directory to store the new namespace equivalence data")
    parser.add_argument(
        "-d",
        required=False,
        action="append",
        help="dataset by prefix; if none specified, all datasets in directory will be run")
    parser.add_argument(
        "-C",
        "--disable_close_matches",
        required=False,
        action="store_true",
        help="disable close matches",
        default=False)
    args = parser.parse_args()
    if os.path.exists(args.n):
        os.chdir(args.n)
    else:
        print('data directory {0} not found!'.format(args.n))

    g = Graph()

    prefix_dict = build_prefix_dict()

    for files in os.listdir("."):
        if files.endswith("parsed_data.pickle"):
            with open(files, 'rb') as f:
                d = pickle.load(f)
                if args.d is None or str(d) in args.d:
                    if isinstance(d, datasets.NamespaceDataSet):
                        make_rdf(d, g, prefix_dict)
                    if isinstance(d, datasets.OrthologyData):
                        get_ortho_matches(d, g, prefix_dict)
                    if isinstance(d, datasets.HistoryDataSet):
                        get_history_data(d, g, prefix_dict)

    if not args.disable_close_matches:
        get_close_matches('BiologicalProcessConcept', g)
        get_close_matches('AbundanceConcept', g)

    print('serializing RDF graph ...')
    g.serialize("testfile.ttl", format='turtle')
