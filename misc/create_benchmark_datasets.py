#!/usr/bin/python2.4
# -*- coding: utf-8 -*-

"""
create_benchmark_datasets.py
This script extracts the datasets used for benchmarking and returns them in JSON and CSV formats for publication on the benchmark captures website.

Created by Shane O'Connor 2015.
Copyright (c) 2015 __UCSF__. All rights reserved.
"""

import sys
import traceback
import datetime
import string
if __name__ == "__main__":
    sys.path.insert(0, "..")
    sys.path.insert(0, "../..")
    sys.path.insert(0, "../updatedb")
import os
import pickle
import tools.colortext as colortext
import tools.deprecated.rosettahelper as rosettahelper
from tools.deprecated.rosettadb import ReusableDatabaseInterface
from tools.fs.fsio import get_file_lines, read_file, write_file
from tools.bio.basics import residue_type_1to3_map, ChainMutation, dssp_elision
from tools.bio.rcsb import retrieve_pdb
from tools.bio.pdb import PDB

from ddgadmin.updatedb.DatasetReferences import Publications, DataSetHomologs

from tools.deprecated.rosettahelper import NUMBER_KJ_IN_KCAL
from ddg.ddglib.ddgdbapi import ddGDatabase, Publicationv2
from ddg.ddglib.ddgobjects import DatasetParser
from ddg.ddglib.db_api import ddG as ddGInterface
import ddg.ddglib.ddgdbapi as ddgdbapi
import ddg.ddglib.ddgobjects as ddgobjects
from tools.biblio.ris import RISEntry
import numpy
import json
import pprint


AAs = set(residue_type_1to3_map.keys())
cached_pdbs = {}
cached_publications = {}
ASEdb = ReusableDatabaseInterface({}, isInnoDB = True, numTries = 1, host = 'guybrush.ucsf.edu', db = 'ASEdb', user = 'oconchus', port = 3306, unix_socket = "/var/lib/mysql/mysql.sock", use_utf = True, passwdfile='asedb.pw')
ddGdb = ddGDatabase(passwd = read_file('ddgdb.pw').strip(), use_utf = False)
ddGdb_utf = ddGDatabase(passwd = read_file('ddgdb.pw').strip(), use_utf = True)
ddG_interface = ddGInterface(passwd = read_file('ddgdb.pw').strip())


JSON_datasets = {
    "AlaScan-GPK_2014/09/25" : dict(
        information = '''
This dataset of experimental alanine scanning data for single mutations was compiled by the Kortemme Lab, UCSF from the following curated datasets:
 - Guerois et al. [1]
 - Potapov et al. [2]
 - Kellogg et al. [3]
and the ProTherm database [4, 5].

The DSSPType values give the DSSP secondary structure type assigned by mkdssp 2.2.1 [6, 7]. The DSSPSimpleType values values
give a simple elision of the DSSP types: G, H, and I become H(elix); E becomes S(heet); and T, B, S, and C become O(ther).
    ''',
        references = {
                1 : 'PMID:12079393',
                2 : 'PMID:19561092',
                3 : 'PMID:21287615',
                4 : 'PMID:9847203',
                5 : 'PMID:16381846',
                6 : 'PMID:21071423',
                7 : 'PMID:6667333',
            }
    ),
    "CuratedProTherm_2014/12/04" : dict(
        information = '''
This dataset of single mutations was compiled by the Kortemme Lab, UCSF from the ProTherm database [1, 2].

The DSSPType column lists the DSSP secondary structure type assigned by mkdssp 2.2.1 [3, 4]. The DSSPSimpleSSType column
lists a simple elision of the DSSP types: G, H, and I become H(elix); E becomes S(heet); and T, B, S, and C become O(ther).
    ''',
        references = {
                1 : 'PMID:9847203',
                2 : 'PMID:16381846',
                3 : 'PMID:21071423',
                4 : 'PMID:6667333',
            }
    ),
    "Guerois_10.1016/S0022-2836(02)00442-4_2002/07/05" : dict(
        information = '''
This dataset of mutations originates from the paper by Guerois et al. [1]. We have made some changes from the original
dataset. In particular, we have tried to match the DDG values to the original publications.

The DSSPType column lists the DSSP secondary structure type assigned by mkdssp 2.2.1 [2, 3]. The DSSPSimpleSSType column
lists a simple elision of the DSSP types: G, H, and I become H(elix); E becomes S(heet); and T, B, S, and C become O(ther).
    ''',
        references = {
                1 : 'PMID:12079393',
                2 : 'PMID:21071423',
                3 : 'PMID:6667333',
            }
    ),
    "Kellogg_10.1002/prot.22921_2010/12/03" : dict(
        information = '''
This dataset of mutations originates from the paper by Kellogg et al. [1]. We have made some changes from the original
dataset. In particular, we have tried to match the DDG values to the original publications.

The DSSPType column lists the DSSP secondary structure type assigned by mkdssp 2.2.1 [2, 3]. The DSSPSimpleSSType column
lists a simple elision of the DSSP types: G, H, and I become H(elix); E becomes S(heet); and T, B, S, and C become O(ther).
    ''',
        references = {
                1 : 'PMID:21287615',
                2 : 'PMID:21071423',
                3 : 'PMID:6667333',
            }
    ),
    "Potapov_10.1093/protein/gzp030_2009/09/01" : dict(
        information = '''
This dataset of mutations originates from the paper by Potapov et al. [1]. We have made some changes from the original
dataset. In particular, we have tried to match the DDG values to the original publications.

The DSSPType column lists the DSSP secondary structure type assigned by mkdssp 2.2.1 [2, 3]. The DSSPSimpleSSType column
lists a simple elision of the DSSP types: G, H, and I become H(elix); E becomes S(heet); and T, B, S, and C become O(ther).
    ''',
        references = {
                1 : 'PMID:19561092',
                2 : 'PMID:21071423',
                3 : 'PMID:6667333',
            }
    ),
}

def read_publications():
    records = ddGdb_utf.execute_select('SELECT * FROM Publication')
    for r in records:
        pubmed_id = ddGdb_utf.execute_select('SELECT * FROM PublicationIdentifier WHERE SourceID=%s AND Type="PMID"', parameters=(r['ID'],))
        if pubmed_id:
            pubmed_id = pubmed_id[0]['ID']
        authors = ddGdb_utf.execute_select('SELECT * FROM PublicationAuthor WHERE PublicationID=%s ORDER BY AuthorOrder', parameters=(r['ID'],))
        authorlist = []
        for a in authors:
            authorlist.append(dict(FirstName = a['FirstName'], MiddleNames = a['MiddleNames'], Surname = a['Surname']))
        pub_details = dict(
            Title = r['Title'],
            Publication = r['Publication'],
            Volume = r['Volume'],
            StartPage = r['StartPage'],
            EndPage = r['EndPage'],
            PublicationYear = r['PublicationYear'],
            PublicationDate = r['PublicationDate'],
            DOI = r['DOI'],
            URL = r['URL'],
            PubMedID = pubmed_id,
            Authors = authorlist,
        )
        if pub_details['PublicationDate']:
            pub_details['PublicationDate'] = pub_details['PublicationDate'].strftime('%Y-%m-%d')

        if not pub_details['URL'] and pub_details['DOI']:
            pub_details['URL'] = 'https://dx.doi.org/%s' % pub_details['DOI']
        cached_publications[r['ID']] = pub_details


def generate_JSON_dataset(dataset_ID, pdb_data, pub_data):

    record_data = {}

    #1LRP
    #1LMB

    # 1 JSON object per dataset record
    failure_count = 0
    records = ddGdb.execute_select('SELECT * FROM DataSetDDG WHERE DataSetID=%s', parameters=(dataset_ID,))
    colortext.warning('Starting with %d records.' % (len(records)))
    mutation_count = {1:0, 2:0, 3:0, 4:0, 5:0}
    for r in records:

        mutation_is_reversed = r['MutationIsReversed'] == 1
        d = dict(
            _DataSetDDGID = r['ID'],
            RecordID = r['RecordNumber'],
            AggregateType = r['AggregateType'],
            DDG = r['PublishedValue'],
            PDBFileID = r['PDBFileID'],
            DerivedMutation = mutation_is_reversed,
        )

        # Parse PDB
        if not(cached_pdbs.get(r['PDBFileID'])):
            cached_pdbs[r['PDBFileID']] = PDB(ddGdb.execute_select('SELECT Content FROM PDBFile WHERE ID=%s', parameters=(r['PDBFileID'],))[0]['Content'])

        # Store PDB data
        PDBResolution = None,
        PDBMethodOfDetermination = None,
        try:
            PDBResolution = cached_pdbs[r['PDBFileID']].get_resolution()
        except: pass
        try:
            PDBMethodOfDetermination = cached_pdbs[r['PDBFileID']].get_techniques()
        except: pass
        pdb_data[r['PDBFileID']] = dict(
            Resolution = PDBResolution,
            MethodOfDetermination = PDBMethodOfDetermination,
        )

        assay_DDGs = ddGdb.execute_select('''
            SELECT *
            FROM DataSetDDGSource
            INNER JOIN ExperimentAssayDDG ON DataSetDDGSource.ExperimentAssayID = ExperimentAssayDDG.ExperimentAssayID AND DataSetDDGSource.Type = ExperimentAssayDDG.Type
            INNER JOIN ExperimentAssay ON ExperimentAssayDDG.ExperimentAssayID = ExperimentAssay.ID
            WHERE DataSetDDGID=%s''', parameters=(r['ID'],))

        ExperimentID = set([a['ExperimentID'] for a in assay_DDGs])
        if len(ExperimentID) != 1:
            colortext.message('%d records passed' % len(record_data))
            # Cases where 1FLV and 1FTG need to be elided
            if sorted(ExperimentID) in ([113699, 113830], [113704, 113832], [113705, 113836]):
                ExperimentID = [sorted(ExperimentID)[0]]
            elif sorted(ExperimentID) in ([112149, 112591],):
                # ExperimentID is used below for mutation details but these agree in this case. 1LZ1, 2BQA
                ExperimentID = [sorted(ExperimentID)[0]]
            elif sorted(ExperimentID) in (
                    [112141, 112583L], [112136, 112578], [112137, 112579], [112142, 112584], [112139, 112581],
                    [112140, 112582], [112146, 112588], [112147, 112589], [112148, 112590]
                ):
                # ExperimentID is used below for mutation details but these agree in this case. 1REX, 2BQA
                ExperimentID = [sorted(ExperimentID)[0]]
            elif sorted(ExperimentID) in ([112227, 112323], [112288, 113039], [111587, 112379]):
                # ExperimentID is used below for mutation details but these agree in this case. 2LZM, 1L63
                ExperimentID = [sorted(ExperimentID)[0]]
            else:
                colortext.warning(
                    '\n'.join(['%(PDBFileID)s %(Chain)s %(WildTypeAA)s %(ResidueID)s %(MutantAA)s' % rii for rii in ddGdb.execute_select('''
                    SELECT * FROM `ExperimentMutation` INNER JOIN Experiment ON Experiment.ID=ExperimentID WHERE `ExperimentID` IN (%s)''' % ','.join(map(str, ExperimentID)))]))
                pprint.pprint(r)
                colortext.error(map(int, ExperimentID))
                #pprint.pprint(assay_DDGs)
                print(sorted(ExperimentID))
        assert(len(ExperimentID) == 1)
        ExperimentID = ExperimentID.pop()
        d['_ExperimentID'] = ExperimentID

        experimental_DDGs = []
        for a in assay_DDGs:
            experimental_DDGs.append(dict(
                DDG = a['Value'],
                DDGType = a['Type'],
                Publication = a['Publication'],
                LocationOfValueInPublication = a['LocationOfValueInPublication'],
                Temperature = a['Temperature'],
                pH= a['pH'],
            ))
            # Store Publication data
            pub_data[a['Publication']] = cached_publications[a['Publication']]
        d['ExperimentalDDGs'] = experimental_DDGs

        # Retrieve mutations
        mutation_records = ddGdb.execute_select('SELECT * FROM ExperimentMutation WHERE ExperimentID=%s ORDER BY ResidueID', parameters=(ExperimentID,))
        if dataset_ID == "AlaScan-GPK_2014/09/25":
            assert(len(mutation_records) == 1)

        mutations = []
        failed_check = False
        mutation_count[len(mutation_records)] += 1
        for mutation in mutation_records:
            mutation_d = {}
            #if ExperimentID == 109911:
            #    d['PDBFileID'] = '1WQ5' # Hack for one 1BKS case

            mutation_d['Chain'] = mutation['Chain']
            mutation_d['ResidueID'] = mutation['ResidueID']
            if mutation_is_reversed:
                mutation_d['MutantAA'] = mutation['WildTypeAA']
                mutation_d['WildTypeAA'] = mutation['MutantAA']
            else:
                mutation_d['WildTypeAA'] = mutation['WildTypeAA']
                mutation_d['MutantAA'] = mutation['MutantAA']

            if dataset_ID == "AlaScan-GPK_2014/09/25":
                if d['PDBFileID'] == '1LMB':
                    mutation_d['Chain'] = '3' # Hack for the PDB replacement 1LRP (3.2A) -> 1LMB (1.8A)
                if d['PDBFileID'] == '1U5P' and int(mutation_d['ResidueID']) < 1600:
                    mutation_d['ResidueID'] = str(int(mutation_d['ResidueID']) + 1762) # Hack for the PDB replacement 1AJ3, NMR -> 1U5P (2A)
            if dataset_ID == "Kellogg_10.1002/prot.22921_2010/12/03":
                if d['PDBFileID'] == '1U5P' and int(mutation_d['ResidueID']) < 1600:
                    mutation_d['ResidueID'] = str(int(mutation_d['ResidueID']) + 1762) # Hack for the PDB replacement 1AJ3, NMR -> 1U5P (2A)

            mutated_residue = ddGdb.execute_select('SELECT * FROM PDBResidue WHERE PDBFileID=%s AND Chain=%s AND ResidueID=%s', parameters=(d['PDBFileID'], mutation_d['Chain'], ResidueID2String(mutation_d['ResidueID'])))
            if len(mutated_residue) == 0:
                colortext.warning('Skipping Experiment #%d (%s) in %s due to missing residue %s.' % (ExperimentID, d['PDBFileID'], dataset_ID, mutation_d['ResidueID']))
                #print('SELECT * FROM PDBResidue WHERE PDBFileID=%s AND Chain=%s AND ResidueID=%s' % (d['PDBFileID'], mutation_d['Chain'], ResidueID2String(mutation_d['ResidueID'])))
                #pprint.pprint(d)
                #pprint.pprint(mutations)
                #pprint.pprint(mutation_d)
                #print(ExperimentID)
                #print(mutated_residue)
                #print(10*'*')
                #print('\n')
                failure_count += 1
                failed_check = True
                break
            assert(len(mutated_residue) == 1)

            mutated_residue = mutated_residue[0]
            mutation_d['DSSPExposure'] = mutated_residue['MonomericExposure']
            mutation_d['DSSPType'] = mutated_residue['MonomericDSSP']
            mutation_d['DSSPSimpleSSType'] = dssp_elision.get(mutation_d['DSSPType'])
            assert(mutation_d['DSSPType'] != None)
            assert(mutation_d['DSSPSimpleSSType'] != None)
            mutations.append(mutation_d)

        if failed_check:
            print('FAILED CHECK')
            continue
        d['Mutations'] = mutations

        if dataset_ID == "Potapov_10.1093/protein/gzp030_2009/09/01":
            key = '%s_%s_%s' % (d['PDBFileID'], '+'.join(['%s:%s:%s' % (mutation_d['Chain'], mutation_d['ResidueID'].strip(), mutation_d['MutantAA']) for mutation_d in mutations]), d['RecordID'])
        else:
            key = '%s_%s' % (d['PDBFileID'], '+'.join(['%s:%s:%s' % (mutation_d['Chain'], mutation_d['ResidueID'].strip(), mutation_d['MutantAA']) for mutation_d in mutations]))

        if record_data.get(key):
            colortext.warning('KEY EXISTS: %s' % key)
            print('Existing record: %s' % pprint.pformat(record_data[key]))
            print('New record: %s' % pprint.pformat(d))
            failure_count += 1
        record_data[key] = d

    colortext.message('Mutation count')
    colortext.warning(pprint.pformat(mutation_count))

    if failure_count > 0:
        colortext.error('Total length of dataset: %d. Failed on %d records.' % (len(record_data), failure_count))
    else:
        colortext.message('Total length of dataset: %d. ' % (len(record_data)))

    record_list = []
    for k, v in sorted(record_data.iteritems()):
        record_list.append(v)

    colortext.message('Adding dataset %s with %d records, %d PDB files, and %d references.' % (dataset_ID, len(record_list), len(pdb_data), len(pub_data)))
    JSON_datasets[dataset_ID]['data'] = record_list



def ChainResidueID2String(chain, residueID):
    '''Takes a chain ID e.g. 'A' and a residueID e.g. '123' or '123A' and returns the 6-character identifier spaced as in the PDB format.'''
    return "%s%s" % (chain, PDB.ResidueID2String(residueID))


def ResidueID2String(residueID):
    '''Takes a residueID e.g. '123' or '123A' and returns the 5-character identifier spaced as in the PDB format.'''
    if residueID.isdigit():
        return "%s " % (residueID.rjust(4))
    else:
        return "%s" % (residueID.rjust(5))


def check_JSON_dataset(dataset_ID):
    # I substitute PDB IDs so this function does a simple check to make sure that the mutations still look okay (this is a simple check - the mutations may not be correct)

    colortext.message('Reading PDB IDs...')
    PDB_ids = set([record['PDBFileID'] for record in JSON_datasets[dataset_ID]['data']])

    colortext.message('Loading %s PDBs...' % len(PDB_ids))
    for PDB_id in PDB_ids:
        if not(cached_pdbs.get(PDB_id)):
            print('Reading %s' % PDB_id)
            colortext.write('.', 'yellow')
            sys.stdout.flush()
            cached_pdbs[PDB_id] = PDB(ddGdb.execute_select('SELECT Content FROM PDBFile WHERE ID=%s', parameters=(PDB_id,))[0]['Content'])
    print('')

    count = 0
    for record in JSON_datasets[dataset_ID]['data']:
        pdb_id = record['PDBFileID']
        p = cached_pdbs[pdb_id]
        #colortext.printf('pdb_id', color='cyan')
        #pprint.pprint(record)
        #pprint.pprint(record['Mutations'])
        for m in record['Mutations']:
            chain_id = m['Chain']
            residue_id = m['ResidueID']
            residue_aa = m['WildTypeAA']
            padded_id = ChainResidueID2String(chain_id, residue_id)
            if p.atom_sequences[chain_id][padded_id].ResidueAA != residue_aa:
                print(pdb_id, chain_id, residue_id, residue_aa)
                print(p.atom_sequences[chain_id][padded_id].ResidueAA, residue_aa)
            assert(p.atom_sequences[chain_id][padded_id].ResidueAA == residue_aa)
        count += 1
    print('Successfully checked %d datapoints.' % count)

def create_dataset_JSON_files():

    todays_date = datetime.date.today().strftime('%Y-%m-%d')
    read_publications()
    pdb_data = {}
    pub_data = {}

    # Add the publications for the datasets
    for k, v in JSON_datasets.iteritems():
        for _, ref in v['references'].iteritems():
            pub_data[ref] = cached_publications[ref]

    #del JSON_datasets["CuratedProTherm_2014/12/04"]
    #del JSON_datasets["Guerois_10.1016/S0022-2836(02)00442-4_2002/07/05"]
    #del JSON_datasets["Potapov_10.1093/protein/gzp030_2009/09/01"]
    #del JSON_datasets["AlaScan-GPK_2014/09/25"]
    #del JSON_datasets["Kellogg_10.1002/prot.22921_2010/12/03"]

    for dataset_ID in JSON_datasets.keys():
        generate_JSON_dataset(dataset_ID, pdb_data, pub_data)
        check_JSON_dataset(dataset_ID)

    max_res = 0
    min_res = 10
    techniques = set()
    for p, v in pdb_data.iteritems():
        techniques.add(v['MethodOfDetermination'])
        if v['Resolution'] != 'N/A':
            max_res = max(max_res, v['Resolution'])
            min_res = min(min_res, v['Resolution'])
    print('Resolutions', min_res, max_res)
    print('Techniques', techniques)

    print(JSON_datasets.keys())
    for k, v in JSON_datasets.iteritems():
        v['version'] = 'This dataset was last updated on %s.' % todays_date
        filename = k.split('_')[0].lower() + '.json'
        x = json.dumps(v, indent=4, sort_keys=True)
        write_file('../rawdata/%s' % filename, x)
    write_file('../rawdata/pdbs.json', json.dumps(pdb_data, indent=4, sort_keys=True))
    write_file('../rawdata/references.json', json.dumps(pub_data, indent=4, sort_keys=True))

def create_dataset_CSV_files():

    csv_lines = []
    csv_lines.append('#' + ','.join(['PDB ID', 'Resolution', 'Techniques']))
    pdbs = json.loads(read_file('../rawdata/pdbs.json'))
    for pdb_id, v in sorted(pdbs.iteritems()):
        csv_lines.append(','.join([pdb_id, str(v['Resolution']), v['MethodOfDetermination']]))
    write_file('../rawdata/pdbs.csv', '\n'.join(csv_lines))

    csv_lines = []
    csv_lines.append('#' + ','.join(['ID', 'Authors', 'Title', 'Publication', 'Volume', 'Issue', 'Date', 'URL']))
    pdbs = json.loads(read_file('../rawdata/references.json'))
    for ref_id, v in sorted(pdbs.iteritems()):
        author_surnames = [a['Surname'] for a in v['Authors']]
        url = ''
        if v['DOI']:
            url = 'https://dx.doi.org/%s' % v['DOI']
        else:
            url = v['URL']
        line = [ref_id, '_'.join(author_surnames), v['Title'], v.get('Publication') or '', v.get('Volume') or '', v.get('Issue') or '', str(v.get('PublicationDate')) or v.get('PublicationYear') or '', url]
        csv_lines.append(','.join(line))
    references_text = '\n'.join(csv_lines)
    import codecs
    f = codecs.open('../rawdata/references.csv', mode="w", encoding="utf-8")
    f.write(references_text)

    for k, v in JSON_datasets.iteritems():
        filename = k.split('_')[0].lower() + '.json'
        json_s = read_file('../rawdata/%s' % filename)
        d = json.loads(json_s)
        csv_lines = []
        print(filename,d.keys())
        for line in d['information'].split('\n'):
            csv_lines.append('# %s' % line)
        for k, v in sorted(d['references'].iteritems()):
            csv_lines.append('# [%s] %s' % (k, v))
        todays_date = datetime.date.today().strftime('%Y-%m-%d')
        csv_lines.append('\n# This dataset was last updated on %s.' % todays_date)
        csv_lines.append('')
        csv_lines.append('# The RecordID below refers to the record ID in the original dataset. When no ID was specified, we added an ID based on the published order of the records.')
        csv_lines.append('# Mutations is an underscore-separated list of mutations. Each mutation takes the form "Chain Wildtype ResidueID Mutant".')
        csv_lines.append('# DDG is the aggregated (mean) DDG value used for analysis.')
        csv_lines.append('# ResidueExposures is an underscore-separated list of exposure values, each one corresponding to its respective mutated position. Each exposure value is based on the solvent accessibility reported by DSSP divided by a maximum solvent accessibility for that residue type and represents whether the residue is buried (0.0) or exposed (1.0).')
        csv_lines.append("# DSSPTypes is an underscore-separated list of DSSP secondary structure assignments, each one corresponding to its respective mutated position.")
        csv_lines.append("# DSSPSimpleTypes is an underscore-separated list of DSSPSimpleType secondary structure assignments, each one corresponding to its respective mutated position.")
        csv_lines.append('# IndividualDDGs lists the individual DDG values which can be used to filter out records with high variance.')
        csv_lines.append('# DerivedMutation is 0 if the record represents an actual set of experiments and 1 if it was derived e.g. if the mutant structure is taken as wildtype and the DDG value is negated. Typically the original records also exist in the dataset so the derived records can introduce a bias.')
        csv_lines.append('# Note: The .json file accompanying this CSV file contains more information, including a list of publications from where the DDG values were taken.')
        csv_lines.append('')

        csv_lines.append('#' + ','.join(['RecordID', 'PDBFileID', 'Mutations', 'DDG', 'ResidueExposures', 'DSSPTypes', 'DSSPSimpleTypes', 'IndividualDDGs', 'DerivedMutation']))
        for record in d['data']:
            mutations = []
            exposures = []
            dssp = []
            ddgs = []
            ddgs_simple = []
            for m in record['Mutations']:
                mutations.append('%(Chain)s %(WildTypeAA)s %(ResidueID)s %(MutantAA)s' % m)
                exposures.append(m['DSSPExposure'])
                dssp.append('%(DSSPType)s' % m)
                ddgs_simple.append('%(DSSPSimpleSSType)s' % m)
            for eddg in record['ExperimentalDDGs']:
                ddgs.append(eddg['DDG'])

            csv_lines.append(','.join([
                str(record['RecordID']),
                record['PDBFileID'],
                '_'.join(mutations),
                str(record['DDG']),
                '_'.join(map(str, exposures)),
                '_'.join(dssp),
                '_'.join(ddgs_simple),
                '_'.join(map(str, ddgs)),
                str(int(record['DerivedMutation'])),
                ]))
        write_file('../rawdata/%s' % filename.replace('.json', '.csv'), '\n'.join(csv_lines))


def dump_pdbs():
    pdbs = json.loads(read_file('../rawdata/pdbs.json'))

    # Sanity check
    for pdb_id, v in sorted(pdbs.iteritems()):
        records = ddGdb.execute_select('SELECT ID FROM PDBFile WHERE ID=%s', parameters=(pdb_id,))
        assert(len(records) == 1)

    # Dump
    for pdb_id, v in sorted(pdbs.iteritems()):
        content = ddGdb.execute_select('SELECT Content FROM PDBFile WHERE ID=%s', parameters=(pdb_id,))[0]['Content']
        write_file('../rawdata/%s.pdb' % pdb_id, content)


def update_public_datasets():
    dsets = ['alascan-gpk.json', 'curatedprotherm.json', 'guerois.json', 'kellogg.json', 'potapov.json']
    source_path = '../rawdata/'
    dest_path = '/home/oconchus/t14benchmarking/ddg/input/json/'
    for dset in dsets:
        assert(os.path.exists(os.path.join(source_path, dset)))
        assert(os.path.exists(os.path.join(dest_path, dset)))
    for dset in dsets:
        print(dset)
        source_set = json.loads(read_file(os.path.join(source_path, dset)))
        dest_set = json.loads(read_file(os.path.join(dest_path, dset)))
        assert(len(source_set['data']) == len(dest_set['data']))
        for x in range(len(source_set['data'])):
            assert(dest_set['data'][x]['RecordID'] == source_set['data'][x]['RecordID'])
            dest_set['data'][x]['DerivedMutation'] = source_set['data'][x]['DerivedMutation']
        write_file(os.path.join(dest_path, dset) + '.new', json.dumps(dest_set, indent=4, sort_keys=True))


if __name__ == '__main__':
    #create_dataset_JSON_files()
    #create_dataset_CSV_files()
    update_public_datasets()
    #dump_pdbs()
    pass
