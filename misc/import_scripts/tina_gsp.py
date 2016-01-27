#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys
import os
import copy
import glob
import datetime
import traceback
import pprint

import pandas

from klab import colortext
from klab.bio.rcsb import retrieve_pdb, download_pdb
from klab.bio.pdb import PDB, coordinate_record_types, chain_record_types
from klab.bio.basics import Residue, Mutation, ChainMutation, generate_all_combinations_of_mutations, SequenceMap
from klab.bio.ligand import LigandMap
from klab.fs.fsio import read_file, write_file, get_file_lines, write_temp_file

sys.path.insert(0, "../../..")

#from ddg.ddglib.ppi_api import get_interface as get_ppi_interface
from ddg.ddglib.ddg_monomer_ppi_api import get_interface as get_ppi_interface
from ddg.ddglib import ddgdbapi, db_api
from ddg.ddglib.import_api import DataImportInterface

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
print(ppi_api.get_unfinished_prediction_ids('ddg_monomer_16_003-zemu-2'))

sys.exit(0)


ppi_api = get_ppi_api()
importer = DataImportInterface.get_interface_with_config_file(cache_dir = '/kortemmelab/data/oconchus/ddgcache')


### Import setup

pdb_file_paths = {}
tina_pdb_objects = {}
rcsb_pdb_objects = {}
tina_pdb_ids = []
tina_pdb_id_to_rcsb_pdb_id = {}
mutations_dataframe = None

def setup():
    global pdb_file_paths  # RCSB PDB_ID -> PDB file
    global rcsb_pdb_objects # RCSB PDB_ID -> PDB object
    global tina_pdb_objects # Tina's PDB_ID -> PDB object
    global tina_pdb_id_to_rcsb_pdb_id # Tina's PDB_ID -> RCSB PDB_ID
    global mutations_dataframe # the updated set of mutations stored as a pandas dataframe with the 'pdb' column converted to uppercase

    # old_mutations_csv is missing some cases but has the mapping from pdb -> partner 1 name, partner 2 name
    old_mutations_csv = os.path.join('temp', 'mutations_Gsp1_old.txt')
    mutations_csv = os.path.join('temp', 'mutations_Gsp1.txt')
    assert(os.path.exists('temp'))
    assert(os.path.exists(old_mutations_csv))
    assert(os.path.exists(mutations_csv))

    df = pandas.read_csv(old_mutations_csv, sep = '\t')
    mutations_dataframe = pandas.read_csv(mutations_csv, sep = '\t')
    mutations_dataframe['pdb'] = mutations_dataframe['pdb'].str.upper()

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

setup()


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

1A2K is in the database with molecules NTF2 (A,B) and RAN/GSP1P (C,D,E). Tina calls this Gsp1|NTF2.
   complex #202 -> 1A2K (C|AB) where Tina uses A|B (which corresponds to C|B after accounting for chain renaming). I will reuse the existing complex and create a new PDBSet.

1K5D is in the database with molecules RAN (A,D,G,J) and RANBP1 (B,E,H,K) and RANGAP (C,F,I,L). Tina calls this Gsp1|RNA1.
   complex #119 -> 1K5D (AB|C). Tina does not use chain B so we will create a new complex.

1I2M is in the database with molecules RAN (A,C) and RCC1 (B,D). Tina calls this Gsp1|SRM1.
   complex #176 -> 1I2M (A|B) where Tina uses A|B (same chain labels). This is an exact match so we will create a new PDBSet.

'''



# ==
#    Step 1
#    Description: Add PDB files
#    Tables: PDBFile, PDBChain, PDBMolecule, PDBMoleculeChain, PDBResidue, PDBLigand, PDBIon
#    Use: import_api.py:DataImportInterface.add_designed_pdb()
# ==
#

user_dataset_name = 'Tina Perica - GSP1 complexes'

def create_mapping_string():
    for pdb_id, tp in sorted(tina_pdb_objects.iteritems()):
        op = rcsb_pdb_objects[tina_pdb_id_to_rcsb_pdb_id[pdb_id]]
        print("'{0}' : dict(".format(pdb_id))
        for c in sorted(tp.atom_sequences.keys()):
            print('\t{0} = ?,'.format(c))
            print(str(tp.atom_sequences[c]))
            print('')
        print('),')


complex_definitions = {
    '1A2K' : dict(
        filepath = 'pdbs/1A2K.pdb', # not used here but this would be the data usually required by e.g. a web API
        rcsb_id = '1A2K',
        db_id = '1A2K_TP0',
        description = 'GSP1 complex from Tina Perica. PDB_REDO was unavailable for this file. Chain A corresponds to chain C in the original RCSB file.',
        params_files = {'G09' : 'temp/pdbs/1A2K.params'},
        LChains = ['A'],
        RChains = ['B'],
        ComplexID = 202,
        chain_mapping = dict(
            A = 'C', # choice of C, D, or E
            B = 'B', # looking at B-factors
        ),
        ligand_mapping = LigandMap.from_tuples_dict({ # Tina's HET code, residue ID -> HET code, RCSB residue ID
            ('G09', 'X   1 ') : ('GDP', 'C 220 '), # PDB columns [17:20], [21:27]
        }),
        ion_mapping = LigandMap.from_tuples_dict({
            ('MG ', 'A 204 ') : ('MG ', 'C 221 '),
        }),
        techniques = "Manual edit",
    ),
}




#tina_pdb_objects : pdb_id -> p
#tina_pdb_id_to_rcsb_pdb_id: pdb_id -> pdb_id
#rcsb_pdb_objects: pdb_id -> p
#tina_to_rcsb_chain_mapping: pdb_id -> chain -> chain

def import_structures():
    ppi_api = get_ppi_api()
    for tina_pdb_id, details in sorted(complex_definitions.iteritems()):
        tina_pdb_id = tina_pdb_id.upper()
        rcsb_pdb_id = tina_pdb_id_to_rcsb_pdb_id[tina_pdb_id]
        assert(details['rcsb_id'] == rcsb_pdb_id)

        tina_pdb_object = tina_pdb_objects[tina_pdb_id]
        rcsb_pdb_object = rcsb_pdb_objects[rcsb_pdb_id]
        tina_db_id = details['db_id']

        assert((tina_db_id != tina_pdb_id) and (tina_db_id != rcsb_pdb_id) and (len(tina_db_id) > 7) and (tina_db_id[4:7] == '_TP'))

        for k, v in details.get('params_files', {}).iteritems():
            details['params_files'][k] = os.path.abspath(v)
        colortext.message('Importing {0} as {1}'.format(tina_pdb_id, tina_db_id))
        importer.add_designed_pdb(tina_pdb_object, tina_db_id, rcsb_pdb_id,
                                  'Tina Perica', details['description'] , 'tina',
                                  chain_mapping = details['chain_mapping'], ligand_mapping = details['ligand_mapping'],
                                  ligand_params_files = details.get('params_files', {}), techniques=details['techniques'])



#import_structures()
sys.exit(0)


sys.exit(0)


tina_to_rcsb_chain_mapping = {
    '1I2M' : dict(
        A = 'A', # choice of A or C
        B = 'B', # choice of B or D
    ),
    '1K5D2' : dict(
        A = 'A', # choice of A, D, G, J
        C = 'C', # choice of C, F, I, L
    ),
    '1QBK' : dict(
        A = 'C',
        B = 'B',
    ),
    '1WA51' : dict(
        A = 'A',
        B = 'B',
    ),
    '1WA52' : dict(
        A = 'A',
        C = 'C',
    ),
    '2BKU' : dict(
        A = 'A', # choice of A, C
        B = 'B', # choice of B, D
    ),
    '3A6P' : dict(
        A = 'C', # choice of C, H   RAN versus the world!
        C = 'A', # choice of A, F
        D = 'D', # choice of D, I
        E = 'E', # choice of E, J
    ),
    '3EA5' : dict(
        A = 'A', # choice of A, C
        B = 'B', # choice of B, D
    ),
    '3ICQ' : dict(
        A = 'B', # choice of B, C   RAN versus the world!
        D = 'D', # choice of D, E
        T = 'T', # choice of T, U
    ),
    '3M1I1' : dict(
        A = 'A',
        B = 'B',
    ),
    '3M1I2' : dict(
        A = 'A',
        C = 'C',
    ),
    '3W3Z' : dict(
        A = 'B',
        B = 'A',
    ),
    '3WYF1' : dict(
        A = 'A', # choice of A, D
        B = 'B', # choice of B, E
    ),
    '4OL0' : dict(
        A = 'A',
        B = 'B',
    ),
}


def add_headers():
    for pdb_id, details in year_2_cases.iteritems():
        source_filepath = details['filepath']
        target_filepath = details['filepath'].replace('.pdb', '.yeast.pdb')
        new_filepath = details['filepath'].replace('.pdb', '.yeast.headers.pdb')
        if read_file(target_filepath).strip().startswith('ATOM'):
            new_content = PDB.replace_headers(read_file(source_filepath), read_file(target_filepath))
            write_file(new_filepath, new_content)





# ==
#    Step 2
#    Description: Add Complex records
#    Tables: PPComplex, PPIPDBSet, PPIPDBPartnerChain
# ==



def create_year_2_complex_records():
    '''todo: For speed, I am doing this manually. This should be folded into an API function.'''

    wt_complex = dict(
        LName = 'Ubiquitin (yeast)',
        LShortName = 'Ubiquitin',
        LHTMLName = 'Ubiquitin (yeast)',
        RName = 'Ubiquitin (yeast)',
        RShortName = 'Ubiquitin',
        RHTMLName = 'Ubiquitin (yeast)',
        FunctionalClassID = 'OX',
        PPDBMFunctionalClassID = 'O',
        PPDBMDifficulty = None,
        IsWildType = True,
        WildTypeComplexID = None,
        Notes = None,
        Warnings = None,
    )
    DDGdb.insertDictIfNew('PPComplex', wt_complex, ['LName', 'RName'])
    wt_complex_id = DDGdb.execute_select('SELECT ID FROM PPComplex WHERE LName=%s and RName=%s', parameters=(wt_complex['LName'], wt_complex['RName']))
    assert(len(wt_complex_id) == 1)
    wt_complex_id = wt_complex_id[0]['ID']

    mut_complex_3H7P = dict(
        LName = 'Ubiquitin (yeast) K63R',
        LShortName = 'Ubiquitin K63R',
        LHTMLName = 'Ubiquitin (yeast) K63R',
        RName = 'Ubiquitin (yeast)',
        RShortName = 'Ubiquitin',
        RHTMLName = 'Ubiquitin (yeast)',
        FunctionalClassID = 'OX',
        PPDBMFunctionalClassID = 'O',
        PPDBMDifficulty = None,
        IsWildType = False,
        WildTypeComplexID = wt_complex_id,
        Notes = None,
        Warnings = None,
    )
    DDGdb.insertDictIfNew('PPComplex', mut_complex_3H7P, ['LName', 'RName'])
    mut_complex_id = DDGdb.execute_select('SELECT ID FROM PPComplex WHERE LName=%s and RName=%s', parameters=(mut_complex_3H7P['LName'], mut_complex_3H7P['RName']))
    assert(len(mut_complex_id) == 1)
    mut_complex_id = mut_complex_id[0]['ID']

    wt_set_number = 0
    for pdb_id, details in sorted(year_2_cases.iteritems()):

        db_pdb_id = 'y' + pdb_id
        if pdb_id == '2W9N':
            db_pdb_id = '2y' + pdb_id

        existing_records = DDGdb.execute_select('''
            SELECT ID
            FROM PPIPDBSet
            INNER JOIN PPIPDBPartnerChain
            WHERE PPIPDBSet.PPComplexID = PPIPDBPartnerChain.PPComplexID
            AND PPIPDBSet.SetNumber = PPIPDBPartnerChain.SetNumber
            AND PPIPDBPartnerChain.PDBFileID=%s''', parameters=(db_pdb_id,))
        if existing_records:
            continue

        complex_id = wt_complex_id
        set_number = None
        if pdb_id == '3H7P':
            complex_id = mut_complex_id
            set_number = 0
        else:
            set_number = wt_set_number
            wt_set_number += 1

        ppi_pdb_set = dict(
            PPComplexID = complex_id,
            SetNumber = set_number,
            IsComplex = True,
            Notes = None,
        )
        ppi_pdb_set_lpartner = dict(
            PPComplexID = complex_id,
            SetNumber = set_number,
            Side = 'L',
            ChainIndex = 0,
            PDBFileID = db_pdb_id,
            Chain = 'A',
            NMRModel = None,
        )
        ppi_pdb_set_rpartner = dict(
            PPComplexID = complex_id,
            SetNumber = set_number,
            Side = 'R',
            ChainIndex = 0,
            PDBFileID = db_pdb_id,
            Chain = 'B',
            NMRModel = None,
        )

        DDGdb.insertDictIfNew('PPIPDBSet', ppi_pdb_set, ['PPComplexID', 'SetNumber'])
        DDGdb.insertDictIfNew('PPIPDBPartnerChain', ppi_pdb_set_lpartner, ['PPComplexID', 'SetNumber', 'Side', 'ChainIndex'])
        DDGdb.insertDictIfNew('PPIPDBPartnerChain', ppi_pdb_set_rpartner, ['PPComplexID', 'SetNumber', 'Side', 'ChainIndex'])


# ==
#    Step 3
#    Description: Add Mutagenesis records
#    Tables: PPMutagenesis, PPMutagenesisMutation, PPMutagenesisPDBMutation
# ==


# ==
#    Step 4
#    Description: Add UserDataSet records
#    Tables: UserDataSet, UserPPDataSetExperiment
# ==


# ==
#    Step 5
#    Description: Add PredictionSet
#    Tables: PredictionSet
#    Use: ppi_api.add_prediction_set()
# ==

prediction_set = 'GSP1 complexes'

# ==
#    Step 6
#    Description: Add Predictions
#    Tables: PredictionPPI
#    Use: ppi_api.add_prediction_run()
# ==


# ==
#    Step 7
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
