#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys
import os
import copy
import glob
import datetime
import traceback
import pprint
import json


import pandas

from klab import colortext
from klab.bio.rcsb import retrieve_pdb, download_pdb
from klab.bio.pdb import PDB, coordinate_record_types, chain_record_types
from klab.bio.basics import Residue, Mutation, ChainMutation, generate_all_combinations_of_mutations, SequenceMap
from klab.bio.ligand import LigandMap
from klab.fs.fsio import read_file, write_file, get_file_lines, write_temp_file

sys.path.insert(0, "../../..")

#from ddg.ddglib.ddg_monomer_ppi_api import get_interface as get_ppi_interface
#from ddg.ddglib.ppi_api import get_interface as get_ppi_interface
from ddg.ddglib.ppi_api import get_interface as get_ppi_interface
from ddg.ddglib import ddgdbapi, db_api
from ddg.ddglib import db_schema as dbmodel
from ddg.ddglib.import_api import DataImportInterface
from klab.db.sqlalchemy_interface import row_to_dict, get_or_create_in_transaction

# This script requires the files in /kortemmelab/data/oconchus/lab_projects/tina/gsp1_ddg/2015November9_ddg.zip
# From this directory:
#     mkdir temp
#     cd temp
#     cp /kortemmelab/data/oconchus/lab_projects/tina/gsp1_ddg/2015November9_ddg.zip .
#     unzip 2015November9_ddg.zip


project_id = 'GSP1'
ppi_api = None

def get_ppi_api():
    global ppi_api
    if not ppi_api:
        ppi_api = get_ppi_interface(read_file('pw'))
    return ppi_api

ppi_api = get_ppi_api()

#ppi_api.merge_prediction_run('ddg_monomer_16-zemu-betanov15', 'ddg_monomer_16-zemu-betanov15v2', create_if_does_not_exist = True)
#ppi_api.merge_prediction_run('ddg_monomer_16-zemu-betanov15', 'ddg_monomer_16-zemu-betanov15v2', create_if_does_not_exist = True)

#ppi_api._ddg_interface._add_file_content('1234567890', db_cursor = None, rm_trailing_line_whitespace = False, forced_mime_type = None)


#print(ppi_api.get_unfinished_prediction_ids('ddg_monomer_16_003-zemu-2'))
#pids = ppi_api.get_prediction_ids_with_scores('ddg_monomer_16_003-zemu-2')
#print(pids)
#print(len(pids), min(pids), max(pids))


### Find all prediction_ids with missing scores and setup rescoring or rescore on the fly
#prediction_ids = ppi_api.get_prediction_ids_with_scores('ddg_monomer_16_003-zemu-2')
#print '%d prediction_ids are yet to be scored' % len(prediction_ids)
#import time
#t1 = time.time()

#ppi_api.get_score_method_details()
#pprint.pprint(ppi_api._ddg_interface.cached_score_method_details)



#ppi_api.get_score_method_details(score_method_id)

#sys.exit(0)

#for prediction_id in list(prediction_ids)[:100]:
#    job_details = ppi_api.get_job_details(prediction_id, truncate_content = 30)
#    pprint.pprint(job_details)
#    break
#    print('.')
#print('{0}s.'.format(time.time() - t1))

# 100 = 7.711s
# 100 new = 8.23s, 8.23s, 7.58s

#sys.exit(0)


#ppi_api = get_ppi_api()
importer = DataImportInterface.get_interface_with_config_file(cache_dir = '/kortemmelab/data/oconchus/ddgcache')


### Import setup

pdb_file_paths = {}
tina_pdb_objects = {}
rcsb_pdb_objects = {}
tina_pdb_ids = []
tina_pdb_id_to_rcsb_pdb_id = {}
mutations_dataframe = None


def setup_mutations_dataframe():
    global mutations_dataframe # the updated set of mutations stored as a pandas dataframe with the 'pdb' column converted to uppercase

    mutations_csv = os.path.join('temp', 'mutations_Gsp1.txt')
    assert(os.path.exists(mutations_csv))
    mutations_dataframe = pandas.read_csv(mutations_csv, sep = '\t')
    mutations_dataframe['pdb'] = mutations_dataframe['pdb'].str.upper()
    mutations_dataframe['DatasetID'] = range(1, len(mutations_dataframe) + 1)
    mutations_dataframe = mutations_dataframe.set_index(['DatasetID'])


def setup():
    global pdb_file_paths  # RCSB PDB_ID -> PDB file
    global rcsb_pdb_objects # RCSB PDB_ID -> PDB object
    global tina_pdb_objects # Tina's PDB_ID -> PDB object
    global tina_pdb_id_to_rcsb_pdb_id # Tina's PDB_ID -> RCSB PDB_ID
    global mutations_dataframe

    if not mutations_dataframe:
        setup_mutations_dataframe()

    # old_mutations_csv is missing some cases but has the mapping from pdb -> partner 1 name, partner 2 name
    old_mutations_csv = os.path.join('temp', 'mutations_Gsp1_old.txt')
    assert(os.path.exists('temp'))
    assert(os.path.exists(old_mutations_csv))

    df = pandas.read_csv(old_mutations_csv, sep = '\t')

    tina_pdb_ids = sorted(set([p for p in df['pdb'].values]))
    rcsb_pdb_ids = set()
    for pdb_id in tina_pdb_ids:
        rcsb_pdb_ids.add(pdb_id[:4])
        tina_pdb_id_to_rcsb_pdb_id[pdb_id] = pdb_id[:4]
    rcsb_pdb_ids = sorted(rcsb_pdb_ids)

    assert(rcsb_pdb_ids == sorted(set([p[:4] for p in mutations_dataframe['pdb'].values])))
    rcsb_file_dir = '../../rawdata'

    for pdb_id in tina_pdb_ids:
        tina_pdb_objects[pdb_id] = PDB.from_filepath(os.path.join('temp', 'pdbs', '{0}.pdb'.format(pdb_id)), parse_ligands = True)

    for pdb_id in rcsb_pdb_ids:
        filename = '{0}.pdb'.format(pdb_id.upper())
        pdb_file_paths[pdb_id.upper()] = os.path.join(rcsb_file_dir, filename)
        pdb_contents = download_pdb(pdb_id, rcsb_file_dir, silent = True, filename = filename)
        p = PDB(pdb_contents, parse_ligands = True)
        rcsb_pdb_objects[pdb_id] = p

    print('\nRosetta files  ({0}) : {1}'.format(str(len(tina_pdb_ids)).rjust(2), ', '.join([s.rjust(5) for s in tina_pdb_ids])))
    print('Original files ({0}) : {1}\n'.format(str(len(rcsb_pdb_ids)).rjust(2), ', '.join([s.rjust(5) for s in rcsb_pdb_ids])))

    ppi_api = get_ppi_api()
    for pdb_id, pdb_file_path in pdb_file_paths.iteritems():
        existing_records = ppi_api.DDG_db.execute_select('SELECT * FROM PDBFile WHERE ID=%s', parameters=(pdb_id,))
        if existing_records:
            colortext.warning('The PDB file {0} exists in the database.'.format(pdb_id))
        complex_ids = ppi_api.search_complexes_by_pdb_id(pdb_id)

        if complex_ids:
            colortext.warning('The PDB file {0} has associated complexes: {1}'.format(pdb_id, ', '.join(map(str, complex_ids))))
    print('')


### Check against existing data
#   Spoiler: There is no experimental binding affinity data at present
###


def print_existing_experimental_data():
    # These PDB files existed in the database before the import so I am interested to see whether any of the experimental
    # data matches the requested predictions
    print('')
    ppi_api = get_ppi_api()
    for pdb_id in ['1A2K', '1K5D', '1I2M']:
        colortext.message(pdb_id)
        complex_ids = ppi_api.search_complexes_by_pdb_id(pdb_id)
        if complex_ids:
            assert(len(complex_ids) == 1)
            complex_id = complex_ids[0]
            colortext.warning('Complex #{0}'.format(complex_id))
            pprint.pprint(ppi_api.get_complex_details(complex_id))

        mutation_records = mutations_dataframe[mutations_dataframe['pdb'].str.contains(pdb_id)]# mutations_dataframe.loc[mutations_dataframe['pdb'][0:4] == pdb_id]
        with pandas.option_context('display.max_rows', None, 'display.max_columns', None):
            print mutation_records

    # There is no experimental binding affinity data at present
    assert(not(ppi_api.DDG_db.execute_select('SELECT * FROM PPMutagenesisPDBMutation WHERE PPComplexID IN (202, 119, 176) ORDER BY PPComplexID, Chain, ResidueID, MutantAA')))


def check_existing_complexes_by_name():
    '''Check whether any of the complexes exist in the database.'''

    # Ran is short for "RAs-related Nuclear protein" and is also known as "GTP-binding nuclear protein Ran"
    ppi_api = get_ppi_api()
    ids = ppi_api.get_complex_ids_matching_protein_name('gsp')
    ids.extend(ppi_api.get_complex_ids_matching_protein_name('ran'))
    ids.extend(ppi_api.get_complex_ids_matching_protein_name('ras'))

    # This gives us these complexes, amongst others:
    #
    # 77
    # Ran GTPase-GDP, Ran GTPase-GDP, Ran GTPase-GDP
    # Importin beta-1 subunit, Importin β1, Importin &beta;1
    #
    # 119
    # Ran GTPase, Ran GTPase, Ran GTPase
    # Ran GAP, Ran GAP, Ran GAP
    #
    # 176
    # Ran GTPase-GDP, Ran GTPase-GDP, Ran GTPase-GDP
    # Regulator of chromosome condensation, RCC1, RCC1
    #
    # 202
    # Ran GTPase-GDP, Ran GTPase-GDP, Ran GTPase-GDP
    # Nuclear transport factor 2, NTF2, NTF2
    #
    # 29
    # Ras GTPase.GDP, Ras GTPase.GDP, Ras GTPase.GDP
    # Ras GAP, Ras GAP, Ras GAP
    #
    # 65
    # Ras GTPase.GTP, H-Ras, H-Ras
    # Son of sevenless-1, Sos, Sos
    #
    # 201
    # Ras GTPase, Ras GTPase, Ras GTPase
    # Phosphoinositide 3-kinase, PI3K, PI3K
    #
    # 280
    # Ras.GNP, Ras.GNP, Ras.GNP
    # RalGDS Ras-interacting domain, RalGDS RID, RalGDS RID

    ids = []
    ids.extend(ppi_api.get_complex_ids_matching_protein_name('importin'))
    ids.extend(ppi_api.get_complex_ids_matching_protein_name('KARYOPHERIN'))
    ids.extend(ppi_api.get_complex_ids_matching_protein_name('TRANSPORTIN'))
    ids.extend(ppi_api.get_complex_ids_matching_protein_name('NTF2'))
    ids.extend(ppi_api.get_complex_ids_matching_protein_name('YRB1P'))
    ids.extend(ppi_api.get_complex_ids_matching_protein_name('RANBP1'))
    ids.extend(ppi_api.get_complex_ids_matching_protein_name('EXP5'))
    ids.extend(ppi_api.get_complex_ids_matching_protein_name('CSE1'))
    ids.extend(ppi_api.get_complex_ids_matching_protein_name('RANGAP'))
    ids.extend(ppi_api.get_complex_ids_matching_protein_name('RANBP2'))
    ids.extend(ppi_api.get_complex_ids_matching_protein_name('RCC1'))

    for id in ids:
        d = ppi_api.get_complex_details(id)
        colortext.warning(id)
        print('{0}, {1}, {2}'.format(d['LName'].encode('utf-8').strip(), d['LShortName'].encode('utf-8').strip(), d['LHTMLName'].encode('utf-8').strip()))
        print('{0}, {1}, {2}'.format(d['RName'].encode('utf-8').strip(), d['RShortName'].encode('utf-8').strip(), d['RHTMLName'].encode('utf-8').strip()))

    # This gives us these complexes:
    #
    # 77
    # Ran GTPase-GDP, Ran GTPase-GDP, Ran GTPase-GDP
    # Importin beta-1 subunit, Importin β1, Importin &beta;1
    #
    # 202
    # Ran GTPase-GDP, Ran GTPase-GDP, Ran GTPase-GDP
    # Nuclear transport factor 2, NTF2, NTF2
    #
    # 176
    # Ran GTPase-GDP, Ran GTPase-GDP, Ran GTPase-GDP
    # Regulator of chromosome condensation, RCC1, RCC1
    #
    # SELECT DISTINCT `PDBFileID` FROM `PPIPDBPartnerChain` WHERE `PPComplexID` IN (77, 202, 176)
    # returns
    # 1F59, 1IBR, 1QG4, 1A12, 1OUN and 1I2M, 1A2K
    #
    # Some of these are unbound. Get the complexes:
    #
    # SELECT DISTINCT `PDBFileID` FROM `PPIPDBPartnerChain`
    # INNER JOIN PPIPDBSet ON PPIPDBPartnerChain.PPComplexID=PPIPDBSet.PPComplexID AND PPIPDBPartnerChain.SetNumber=PPIPDBSet.SetNumber
    # WHERE PPIPDBPartnerChain.PPComplexID IN (77, 202, 176) AND IsComplex=1
    #
    # returns only three hits:
    #  complex #77  -> 1IBR (A|B);
    #  complex #176 -> 1I2M (A|B) where Tina uses A|B (chains may be renamed); and
    #  complex #202 -> 1A2K (C|AB) where Tina uses A|B (chains may be renamed).
    #
    # We also have:
    #  complex #119 -> 1K5D (AB|C) where Tina uses A|B
    #
    # 1IBR -> Ran (human)|Importin β1 (human)
    # Tina has:
    #    2BKU -> RAN (dog)|Importin β1 (yeast)
    #    3EA5 -> RAN (human)|Importin β1 (yeast)
    # 3EA5 and 1IBR do not match on chains B at all and have one mutation in chain A
    # Similarly for 2BKU and 1IBR.
    #
    # However what came out of this is that 3EA5 and 2BKU are related i.e. that RAN is almost the same sequence in both.
    # The only difference is one mutation in chain A: index 40, A->P and that 3EA5 has a longer sequence for chain A
    #

    colortext.message('\n\n1IBR')
    p1 = PDB(retrieve_pdb('1IBR'))
    pprint.pprint(p1.seqres_sequences)
    colortext.message('\n\n2BKU')
    p2 = PDB(retrieve_pdb('2BKU'))
    pprint.pprint(p2.seqres_sequences)
    a1 = str(p1.seqres_sequences['A'])
    a2 = str(p2.seqres_sequences['A'])

    #3EA5
    a1 = 'MAAQGEPQVQFKLVLVGDGGTGKTTFVKRHLTGEFEKKYVATLGVEVHPLVFHTNRGPIKFNVWDTAGQEKFGGLRDGYYIQAQCAIIMFDVTSRVTYKNVPNWHRDLVRVCENIPIVLCGNKVDIKDRKVKAKSIVFHRKKNLQYYDISAKSNYNFEKPFLWLARKLIGDPNLEFVAMPCLAPPEVVMDPALAAQYEHDLEVAQTTALPDEDDDL'
    a1 = 'MSTAEFAQLLENSILSPDQNIRLTSETQLKKLSNDNFLQFAGLSSQVLIDENTKLEGRILAALTLKNELVSKDSVKTQQFAQRWITQVSPEAKNQIKTNALTALVSIEPRIANAAAQLIAAIADIELPHGAWPELMKIMVDNTGAEQPENVKRASLLALGYMCESADPQSQALVSSSNNILIAIVQGAQSTETSKAVRLAALNALADSLIFIKNNMEREGERNYLMQVVCEATQAEDIEVQAAAFGCLCKIMSKYYTFMKPYMEQALYALTIATMKSPNDKVASMTVEFWSTICEEEIDIAYELAQFPQSPLQSYNFALSSIKDVVPNLLNLLTRQNEDPEDDDWNVSMSAGACLQLFAQNCGNHILEPVLEFVEQNITADNWRNREAAVMAFGSIMDGPDKVQRTYYVHQALPSILNLMNDQSLQVKETTAWCIGRIADSVAESIDPQQHLPGVVQACLIGLQDHPKVATNCSWTIINLVEQLAEATPSPIYNFYPALVDGLIGAANRIDNEFNARASAFSALTTMVEYATDTVAETSASISTFVMDKLGQTMSVDENQLTLEDAQSLQELQSNILTVLAAVIRKSPSSVEPVADMLMGLFFRLLEKKDSAFIEDDVFYAISALAASLGKGFEKYLETFSPYLLKALNQVDSPVSITAVGFIADISNSLEEDFRRYSDAMMNVLAQMISNPNARRELKPAVLSVFGDIASNIGADFIPYLNDIMALCVAAQNTKPENGTLEALDYQIKVLEAVLDAYVGIVAGLHDKPEALFPYVGTIFQFIAQVAEDPQLYSEDATSRAAVGLIGDIAAMFPDGSIKQFYGQDWVIDYIKRTRSGQLFSQATKDTARWAREQQKRQLSL'
    #2BKU
    a2 = 'MAAQGEPQVQFKLVLVGDGGTGKTTFVKRHLTGEFEKKYVPTLGVEVHPLVFHTNRGPIKFNVWDTAGQEKFGGLRDGYYIQAQCAIIMFDVTSRVTYKNVPNWHRDLVRVCENIPIVLCGNKVDIKDRKVKAKSIVFHRKKNLQYYDISAKSNYNFEKPFLWLARKLIGDPNLEFV'
    a2 = 'MSTAEFAQLLENSILSPDQNIRLTSETQLKKLSNDNFLQFAGLSSQVLIDENTKLEGRILAALTLKNELVSKDSVKTQQFAQRWITQVSPEAKNQIKTNALTALVSIEPRIANAAAQLIAAIADIELPHGAWPELMKIMVDNTGAEQPENVKRASLLALGYMCESADPQSQALVSSSNNILIAIVQGAQSTETSKAVRLAALNALADSLIFIKNNMEREGERNYLMQVVCEATQAEDIEVQAAAFGCLCKIMSKYYTFMKPYMEQALYALTIATMKSPNDKVASMTVEFWSTICEEEIDIAYELAQFPQSPLQSYNFALSSIKDVVPNLLNLLTRQNEDPEDDDWNVSMSAGACLQLFAQNCGNHILEPVLEFVEQNITADNWRNREAAVMAFGSIMDGPDKVQRTYYVHQALPSILNLMNDQSLQVKETTAWCIGRIADSVAESIDPQQHLPGVVQACLIGLQDHPKVATNCSWTIINLVEQLAEATPSPIYNFYPALVDGLIGAANRIDNEFNARASAFSALTTMVEYATDTVAETSASISTFVMDKLGQTMSVDENQLTLEDAQSLQELQSNILTVLAAVIRKSPSSVEPVADMLMGLFFRLLEKKDSAFIEDDVFYAISALAASLGKGFEKYLETFSPYLLKALNQVDSPVSITAVGFIADISNSLEEDFRRYSDAMMNVLAQMISNPNARRELKPAVLSVFGDIASNIGADFIPYLNDIMALCVAAQNTKPENGTLEALDYQIKVLEAVLDAYVGIVAGLHDKPEALFPYVGTIFQFIAQVAEDPQLYSEDATSRAAVGLIGDIAAMFPDGSIKQFYGQDWVIDYIKRTRSGQLFSQATKDTARWAREQQKRQLSL'
    print(a1 == a2)
    if not a1 == a2:
        # horribly inefficient (casting to str each time) but not worth rewriting
        assert(len(a1) == len(a2))
        for x in range(len(a1)):
            if str(a1)[x] != str(a2)[x]:
                print(x, str(a1)[x], str(a2)[x])
        # one mutation A->C near the end of the sequence: VAMPALAP -> VAMPCLAP

    assert(str(p1.seqres_sequences['A']) == str(p1.seqres_sequences['C']))
    assert(str(p1.seqres_sequences['B']) == str(p1.seqres_sequences['D']))
    assert(str(p2.seqres_sequences['A']) == str(p2.seqres_sequences['C']))
    assert(str(p2.seqres_sequences['B']) == str(p2.seqres_sequences['D']))
    print('')

'''
Did Tina change the chain letters of any PDB files?
Yes. Tina always has gsp1 as chain A?
The mapping is below.

At this point, we have a mapping from three of Tina's complexes to entries in the database but with no corresponding experimental data.

1I2M is in the database with molecules RAN (A,C) and RCC1 (B,D). Tina calls this Gsp1|SRM1.
   complex #176 -> 1I2M (A|B) where Tina uses A|B (same chain labels). This is an exact match so we will use this complex (but we still need to add a new PDBSet since this is a different PDB file).

1A2K is in the database with molecules NTF2 (A,B) and RAN/GSP1P (C,D,E). Tina calls this Gsp1|NTF2.
   complex #202 -> 1A2K (C|AB) where Tina uses A|B (which corresponds to C|B after accounting for chain renaming). I will reuse the existing complex and create a new PDBSet.

1K5D is in the database with molecules RAN (A,D,G,J) and RANBP1 (B,E,H,K) and RANGAP (C,F,I,L). Tina calls this Gsp1|RNA1.
   complex #119 -> 1K5D (AB|C). Tina does not use chain B so we will create a new complex.

'''



# ==
#    Step 1
#    Description: Add PDB files and complexes
#    Use: import_api.py:DataImportInterface.add_complex_structure_pair() - this calls add_designed_pdb() and add_complex()
# ==
#

user_dataset_name = 'Tina Perica - GSP1 complexes'

def create_mapping_string():
    # todo: use this to create a function to create a skeleton JSON file given a list of complexes
    for pdb_id, tp in sorted(tina_pdb_objects.iteritems()):
        op = rcsb_pdb_objects[tina_pdb_id_to_rcsb_pdb_id[pdb_id]]
        print("'{0}' : dict(".format(pdb_id))
        for c in sorted(tp.atom_sequences.keys()):
            print('\t{0} = ?,'.format(c))
            print(str(tp.atom_sequences[c]))
            print('')
        print('),')


def create_project_pdb_records():
    ppi_api = get_ppi_api()
    complex_definitions = json.loads(read_file('tinas_complexes.json'))
    for k, v in complex_definitions.iteritems():
        ppi_api.associate_pdb_file_with_project(v['Structure']['db_id'], u'GSP1')


def import_structures():
    setup()
    ppi_api = get_ppi_api()
    complex_definitions = json.loads(read_file('tinas_complexes.json'))
    for tina_pdb_id, complex_structure_definition_pair in sorted(complex_definitions.iteritems()):
        #if tina_pdb_id != '1WA52':
        #    continue
        colortext.warning(tina_pdb_id)
        del complex_structure_definition_pair['Structure']['file_path']
        complex_structure_definition_pair['Structure']['pdb_object'] = tina_pdb_objects[tina_pdb_id]
        pdb_set = ppi_api.add_complex_structure_pair(complex_structure_definition_pair, keywords = ['GSP1'],
                                                     force = True, trust_database_content = False, allow_missing_params_files = False, debug = False)
        if pdb_set['success'] == False:
            print(pdb_set['error'])
            if 'possible_matches' in pdb_set:
                for d in pdb_set['possible_matches']:
                    colortext.warning(d['ID'])
                    print('{0}, {1}, {2}'.format(d['LName'].encode('utf-8').strip(), d['LShortName'].encode('utf-8').strip(), d['LHTMLName'].encode('utf-8').strip()))
                    print('{0}, {1}, {2}'.format(d['RName'].encode('utf-8').strip(), d['RShortName'].encode('utf-8').strip(), d['RHTMLName'].encode('utf-8').strip()))

    create_project_pdb_records()


# ==
#    Step 2
#    Description: Add Mutagenesis and UserDataSet records
#    Tables: PPMutagenesis, PPMutagenesisMutation, PPMutagenesisPDBMutation
# ==


def import_mutageneses():
    setup_mutations_dataframe()
    ppi_api = get_ppi_api()
    
    complex_definitions = json.loads(read_file('tinas_complexes.json'))

    # Determine the mapping from PDB ID to complex ID
    pdb_id_to_database_id = {}
    for index, r in mutations_dataframe.iterrows():
        pdb_id = r['pdb']
        db_id = complex_definitions[pdb_id]['Structure']['db_id']
        if pdb_id_to_database_id.get(pdb_id):
            assert(pdb_id_to_database_id[pdb_id] == db_id)
        pdb_id_to_database_id[pdb_id] = db_id

    pdb_id_to_complex_id = {}
    for pdb_id, db_id in sorted(pdb_id_to_database_id.iteritems()):
        results = ppi_api.DDG_db.execute_select('SELECT DISTINCT PPComplexID, SetNumber FROM PPIPDBPartnerChain WHERE PDBFileID=%s', parameters=(db_id,))
        assert(len(results) == 1)
        pdb_id_to_complex_id[pdb_id] = dict(PPComplexID = results[0]['PPComplexID'], SetNumber = results[0]['SetNumber'])

    pdb_residues = {}
    for db_id in pdb_id_to_database_id.values():
        pdb_residues[db_id] = {}
        for r in ppi_api.DDG_db.execute_select('SELECT Chain, ResidueID, ResidueAA FROM PDBResidue WHERE PDBFileID=%s', parameters=(db_id,)):
            pdb_residues[db_id][r['Chain']] = pdb_residues[db_id].get(r['Chain'], {})
            pdb_residues[db_id][r['Chain']][r['ResidueID']] = r['ResidueAA']

    assert(len(pdb_id_to_complex_id) == 15)

    user_data_set_text_id = 'RAN-GSP'
    ppi_api.add_user_dataset('oconchus', user_data_set_text_id, "Tina's dataset for RAN/GSP1 complexes.")

    user_dataset_cases = []
    for index, r in mutations_dataframe.iterrows():
        pdb_id = r['pdb']
        database_pdb_id = pdb_id_to_database_id[pdb_id]
        dataset_id = index
        pdb_id = r['pdb']
        complex_definition = complex_definitions[pdb_id]

        # all the mutations are on chain1 (which is always chain A)
        chain_id = 'A'
        residue_id = str(r['pdb_res_num'])
        wildtype_aa = pdb_residues[database_pdb_id][chain_id][PDB.ResidueID2String(residue_id)]
        mutant_aa = r['mutation']
        assert(wildtype_aa != mutant_aa)

        case_details = dict(

            # These records are used to create a PPMutagenesis record and the associated mutagenesis details

            Mutagenesis = dict(
                RecognizableString = 'TinaGSP_{0}'.format(dataset_id),
                PPComplexID = pdb_id_to_complex_id[pdb_id]['PPComplexID'],
            ),

            Mutations = [
                # There is one dict per mutation
                dict(
                    MutagenesisMutation = dict(
                        # PPMutagenesisID will be filled in when the PPMutagenesis record is created.
                        RecordKey = '{0} {1}{2}{3}'.format(chain_id, wildtype_aa, residue_id.strip(), mutant_aa),
                        ProteinID = None, # todo
                        ResidueIndex = None, # todo
                        WildTypeAA = wildtype_aa,
                        MutantAA = mutant_aa,
                    ),
                    MutagenesisPDBMutation = dict(
                        # PPMutagenesisID and PPMutagenesisMutationID will be filled in when the PPMutagenesisMutation record is created.
                        # PPComplexID is taken from the PPMutagenesis section. WildTypeAA and MutantAA are taken from the PPMutagenesisMutation section.
                        SetNumber = pdb_id_to_complex_id[pdb_id]['SetNumber'],
                        PDBFileID = database_pdb_id,
                        Chain = chain_id,
                        ResidueID = residue_id,
                    ),
                ),
            ],

            # This field is used to create the UserPPDataSetExperiment record. All other fields can be derived from the above.
            # Note: We use the human-readable label here. The database ID is retrieved using e.g. ppi_api.get_defined_user_datasets()[<UserDataSetTextID>]['ID']
            UserDataSetTextID = user_data_set_text_id,
        )
        user_dataset_cases.append(case_details)

    colortext.porange('Creating the UserDataSet cases')
    user_dataset_name_to_id_map = {}
    tsession = ppi_api.get_session(new_session = True)
    try:
        for user_dataset_case in user_dataset_cases:
            ppi_api.add_user_dataset_case(tsession, user_dataset_case, user_dataset_name_to_id_map = user_dataset_name_to_id_map)

        print('\n\nSuccess')
        tsession.commit()
        #tsession.rollback()
        tsession.close()
    except Exception, e:
        colortext.error('\n\nFailure: An error occurred.')
        colortext.warning(str(e))
        colortext.warning(traceback.format_exc())
        tsession.rollback()
        tsession.close()

# ==
#    Step 3
#    Description: Add PredictionSet and predictions
#    Tables: PredictionSet
#    Use: ppi_api.add_prediction_set()
#    Use: ppi_api.add_prediction_run()
# ==

prediction_set = 'GSP1 complexes'


# ==
#    Step 4
#    Description: Test prediction set
#    Use: ppi_api.get_queued_jobs()
# ==


def test_prediction_set():
    c = 0
    counts = {}
    ppi_api = get_ppi_api()
    for j in ppi_api.get_queued_jobs(prediction_set, order_by = 'Cost', order_order_asc = False, include_files = False, truncate_content = None):
        counts[j['Structure']['PDBFileID']] = counts.get(j['Structure']['PDBFileID'], 0)
        counts[j['Structure']['PDBFileID']] += 1
        c += 1
    colortext.warning('Counts by PDB ID:')
    pprint.pprint(counts)
    colortext.warning('Total count: {0}'.format(c))
