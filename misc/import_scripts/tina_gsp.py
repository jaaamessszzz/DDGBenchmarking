#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys
import os
import copy
import glob
import datetime
import traceback
import pprint
#sys.path.insert(0, "../../../klab")
#sys.path.insert(0, "../../..")

import pandas

from klab import colortext
from klab.bio.rcsb import retrieve_pdb, download_pdb
from klab.bio.pdb import PDB, coordinate_record_types, chain_record_types
from klab.bio.basics import Residue, Mutation, ChainMutation, generate_all_combinations_of_mutations, SequenceMap
from klab.fs.fsio import read_file, write_file, get_file_lines, write_temp_file

from ddg.ddglib.ppi_api import get_interface as get_ppi_interface
from ddg.ddglib import ddgdbapi, db_api

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

pdb_file_paths = {}
tina_pdb_objects = {}
rcsb_pdb_objects = {}
tina_pdb_ids = []
tina_pdb_id_to_rcsb_pdb_id = {}
mutations_dataframe = None

def setup():
    global pdb_file_paths
    global rcsb_pdb_objects
    global tina_pdb_objects
    global tina_pdb_id_to_rcsb_pdb_id
    global mutations_dataframe

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
        tina_pdb_objects[pdb_id] = PDB.from_filepath(os.path.join('temp', 'pdbs', '{0}.pdb'.format(pdb_id)))

    for pdb_id in rcsb_pdb_ids:
        filename = '{0}.pdb'.format(pdb_id.upper())
        pdb_file_paths[pdb_id.upper()] = os.path.join(rcsb_file_dir, filename)
        pdb_contents = download_pdb(pdb_id, rcsb_file_dir, silent = True, filename = filename)
        p = PDB(pdb_contents)
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
At this point, we have a mapping from three of Tina's complexes to entries in the database but with no corresponding experimental data.

1A2K is in the database with molecules NTF2 (A,B) and RAN/GSP1P (C,D,E). Tina calls this Gsp1|NTF2.
   complex #202 -> 1A2K (C|AB) where Tina uses A|B (chains may be renamed).

1K5D is in the database with molecules RAN (A,D,G,J) and RANBP1 (B,E,H,K) and RANGAP (C,F,I,L). Tina calls this Gsp1|RNA1.
   complex #119 -> 1K5D

1I2M is in the database with molecules RAN (A,C) and RCC1 (B,D). Tina calls this Gsp1|SRM1.
   complex #176 -> 1I2M (A|B) where Tina uses A|B (chains may be renamed); and

did Tina change the chain letters of any PDB files?
yes. Tina always has gsp1 as chain A?
'''




### Create the PDB records

##########################
year_2_cases = {
    '3NOB' : dict(
        filepath = '/kortemmelab/data/oconchus/PUBS/year2/K11_3NOB.pdb',
        conserved = [dict(ResidueID = '11', ResidueAA = 'K')],
        LChains = ['A'],
        RChains = ['B'],
        prediction_set = 'Ubiquitin scan: K11_3NOB p16',
    ),
    '4S22' : dict(
        filepath = '/kortemmelab/data/oconchus/PUBS/year2/K29_4S22.pdb',
        conserved = [dict(ResidueID = '29', ResidueAA = 'K')],
        LChains = ['A'],
        RChains = ['B'],
        prediction_set = 'Ubiquitin scan: K29_4S22 p16',
    ),
    '1AAR' : dict(
        filepath = '/kortemmelab/data/oconchus/PUBS/year2/K48_1AAR.pdb',
        conserved = [dict(ResidueID = '48', ResidueAA = 'K')],
        LChains = ['A'],
        RChains = ['B'],
        prediction_set = 'Ubiquitin scan: K48_1AAR p16',
    ),
    '2XK5' : dict(
        filepath = '/kortemmelab/data/oconchus/PUBS/year2/K6_2XK5.pdb',
        conserved = [dict(ResidueID = '6', ResidueAA = 'K')],
        LChains = ['A'],
        RChains = ['B'],
        prediction_set = 'Ubiquitin scan: K6_2XK5 p16',
    ),
    '3H7P' : dict(
        filepath = '/kortemmelab/data/oconchus/PUBS/year2/K63_3H7P.pdb',
        conserved = [dict(ResidueID = '63', ResidueAA = 'K')],
        LChains = ['A'],
        RChains = ['B'],
        prediction_set = 'Ubiquitin scan: K63_3H7P p16',
    ),
    '2W9N' : dict(
        filepath = '/kortemmelab/data/oconchus/PUBS/year2/M1_2W9N.pdb',
        conserved = [dict(ResidueID = '1', ResidueAA = 'M')],
        LChains = ['A'],
        RChains = ['B'],
        prediction_set = 'Ubiquitin scan: M1_2W9N p16',
    ),
}
#####################

def create_mapping_string():
    for pdb_id, tp in sorted(tina_pdb_objects.iteritems()):
        op = rcsb_pdb_objects[tina_pdb_id_to_rcsb_pdb_id[pdb_id]]
        print("'{0}' : dict(".format(pdb_id))
        for c in sorted(tp.atom_sequences.keys()):
            print('\t{0} = ?,'.format(c))
            print(str(tp.atom_sequences[c]))
            print('')
        print('),')

tina_to_rcsb_chain_mapping = {
    '1A2K' : dict(
        A = 'C', # choice of C, D, or E
        B = 'A', # choice of A or B
    ),
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
        A = 'C', # choice of C, H RAN versus the world!
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

#tina_pdb_objects : pdb_id -> p
#tina_pdb_id_to_rcsb_pdb_id: pdb_id -> pdb_id
#rcsb_pdb_objects: pdb_id -> p
#tina_to_rcsb_chain_mapping: pdb_id -> chain -> chain

#1. Update the PDB function to allow renaming of chains i.e. allow it to take a non-RCSB PDB file, a mapping from new chain letters to RCSB letters, and
#use the header information e.g. molecule etc. renamed for the artificial structure



def add_headers():
    for pdb_id, details in year_2_cases.iteritems():
        source_filepath = details['filepath']
        target_filepath = details['filepath'].replace('.pdb', '.yeast.pdb')
        new_filepath = details['filepath'].replace('.pdb', '.yeast.headers.pdb')
        if read_file(target_filepath).strip().startswith('ATOM'):
            new_content = PDB.replace_headers(read_file(source_filepath), read_file(target_filepath))
            write_file(new_filepath, new_content)

# ==
#    Step 1
#    Description: Add PDB files
#    Tables: PDBFile, PDBChain, PDBMolecule, PDBMoleculeChain, PDBResidue
#    Use: ppi_api.add_PDB_to_database()
# ==
#

def import_year_2_structures():
    ppi_api = get_ppi_interface(read_file('pw'))
    for tina_pdb_id, p in sorted(tina_pdb_objects.iteritems()):
        db_name = 'tp{0}'.format(tina_pdb_id)
        colortext.message('Importing {0} as {1}'.format(tina_pdb_id, db_name))

        # The UniProt mapping expects the reference PDB file to exist in the database so we will add it
        rcsb_pdb_id = tina_pdb_id_to_rcsb_pdb_id[tina_pdb_id]
        results = ppi_api.DDG_db.execute_select('SELECT * FROM PDBFile WHERE ID=%s', parameters=(rcsb_pdb_id,))
        if len(results) == 0:
            #fname = write_temp_file('/tmp', retrieve_pdb(rcsb_pdb_id), suffix = '.pdb')
            #print(fname)
            ppi_api.add_RCSB_pdb_to_database(pdb_id = rcsb_pdb_id)
            continue
            sys.exit(0)
            ppi_api.add_PDB_to_database(
                filepath = fname,
                pdbID = rcsb_pdb_id,
                force = True,
                file_source = 'RCSB',
            )
            os.remove(fname)
        continue
        ppi_api.add_PDB_to_database(
            filepath = details['db_filepath'],
            pdbID = 'y' + pdb_id,
            derived_from = pdb_id,
            force = True,
            file_source = 'Rosetta',
            notes = "Created for the PUBS class at UCSF. Contact David Mavor, Kyle Barlow, Samuel Thompson, or Shane O'Connor for more details. This file is derived from an RCSB human ubiquitin structure but has been altered using the Rosetta fixbb application to use the yeast ubiquitin sequence. HETATM records and chains may also have been removed.",
            allow_missing_molecules = True
        )


import_year_2_structures()
sys.exit(0)


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
