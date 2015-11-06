import sys
import os
import copy
import glob
import datetime
import traceback
import pprint
sys.path.insert(0, "../klab")


from klab import colortext
from klab.bio.rcsb import retrieve_pdb
from klab.bio.pdb import PDB, chain_record_types
from klab.bio.basics import Residue, Mutation, ChainMutation, generate_all_combinations_of_mutations
from klab.fs.fsio import read_file, write_file, get_file_lines, write_temp_file

from ddglib.ppi_api import get_interface as get_ppi_interface
from ddglib import ddgdbapi, db_api

if __name__ == '__main__':
    DDGdb = ddgdbapi.ddGDatabase()
    ddG_connection = db_api.ddG()

all_wildtype_mutations = '''
# All mutations are on chain A
A83P, I86L, S87A
L95I, I98V
A83P, I86L, S87A, L95I, I98V
'''


ubiquitin_chains = [
    ('1ubq', 'A', 'Ubiquitin scan: 1UBQ p16'),
    ('uby_1UBQ', 'A', 'Ubiquitin scan: 1UBQ_yeast p16'),
    ('ub_SH3', 'A', 'Ubiquitin scan: SH3 p16'),
    ('uby_SH3', 'A', 'Ubiquitin scan: SH3_yeast p16'),
    ('ub_CUE', 'A', 'Ubiquitin scan: CUE p16'),
    ('uby_CUE', 'A', 'Ubiquitin scan: CUE_yeast p16'),
    ('ub_OTU', 'A', 'Ubiquitin scan: OTU p16'),
    ('uby_OTU', 'A', 'Ubiquitin scan: OTU_yeast p16'),
    ('ub_RPN13', 'A', 'Ubiquitin scan: RPN13 p16'),
    ('uby_RPN13', 'A', 'Ubiquitin scan: RPN13_yeast p16'),
    ('ub_UQcon', 'A', 'Ubiquitin scan: UQ_con p16'),
    ('uby_UQcon', 'A', 'Ubiquitin scan: UQ_con_yeast p16'),
]

human_sequence = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRG'
yeast_sequence = 'MQIFVKTLTGKTITLEVESSDTIDNVKSKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRG'

def get_data():
    for ubiquitin_chain in ubiquitin_chains:
        results = DDGdb.execute_select('''
            SELECT PredictionSet, PDBFileID, COUNT( Prediction.ID ) AS NumRecords
            FROM `Prediction`
            INNER JOIN Experiment ON ExperimentID=Experiment.ID
            WHERE PredictionSet=%s
            GROUP BY PredictionSet
            ''', parameters = (ubiquitin_chain[2],))
        assert(len(results) == 1)
        assert(1425 <= results[0]['NumRecords'] <= 1445)

        cases = []
        results = DDGdb.execute_select('''
            SELECT Prediction.ID AS PredictionID, ProtocolID, Prediction.ExperimentID, PredictionSet, PDBFileID,
            Status, Errors
            FROM `Prediction`
            INNER JOIN Experiment ON ExperimentID=Experiment.ID
            WHERE PredictionSet=%s''', parameters = (ubiquitin_chain[2],))
        native = True
        organism = 'Human'
        if ubiquitin_chain[0].find('uby_') != 0:
            organism = 'Yeast'
            native = False
        case_name = ubiquitin_chain[0].replace('ub_', '').replace('uby_', '')
        for r in results:
            #print(r)
            mutations = DDGdb.execute_select('SELECT * FROM ExperimentMutation WHERE ExperimentID=%s', parameters=(r['ExperimentID'],))
            assert(len(mutations) == 1)
            mutation = mutations[0]
            assert(r['Status'] == 'done' and r['Errors'] == None)
            prediction_id = r['PredictionID']
            archive = '/kortemmelab/shared/DDG/jobs/{0}.zip'.format(prediction_id)
            assert(os.path.exists(archive))
            cases.append(dict(
                case_name = case_name,
                organism = organism,
                native = native,
                prediction_id = r['PredictionID'],
                protocol_id = r['ProtocolID'],
                experiment_id = r['ExperimentID'],
                prediction_set = r['PredictionSet'],
                pdb_file_id = r['PDBFileID'],
                zip_filepath = archive,
                chain = mutation['Chain'],
                mutant_aa = mutation['MutantAA'],
                wildtype_aa = mutation['WildTypeAA'],
                residue_id = PDB.ResidueID2String(mutation['ResidueID']),
            ))

        # run multiprocessing over cases
        # if new score type already stored in PredictionStructureScore, continue
        # else
        #     compute score from zip
        #     add record to PredictionStructureScore using new ScoreMethod ID and appropriate ScoreType & StructureID
        #     for run, StructureID should be -1 and ScoreType should be DDG
        #     for structure/pair, StructureID should be 1-50 depending on filename and ScoreType should be something sensible - at present I only use Mutant and WildType because I treated these runs as monomers


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

def setup_year_2_yeast_structures():

    hetatm_types = set()
    for pdb_id, details in year_2_cases.iteritems():
        print(pdb_id)
        filepath = details['filepath']
        pruned_filepath = os.path.splitext(filepath)[0] + '.pruned.pdb'
        assert(os.path.exists(filepath))
        p = PDB.from_filepath(filepath)

        resfile = []
        if pdb_id != '2W9N':
            assert('A' in p.atom_sequences and 'B' in p.atom_sequences)

            # Write new file
            pruned_lines = []
            for l in get_file_lines(filepath):
                if l[:6] == 'HETATM':
                    hetatm_types.add(l[17:20])
                if l[:6].strip() in chain_record_types and l[21] != 'A' and l[21] != 'B':
                    continue
                if l[:6] != 'HETATM' or (l[17:20] == 'HOH'): # it seems fine to remove ZN ions in these cases
                    pruned_lines.append(l)
            write_file(pruned_filepath, '\n'.join(pruned_lines))
            p = PDB.from_filepath(pruned_filepath)
            assert(sorted(p.atom_sequences.keys()) == ['A', 'B'])

            chain_a_seq_o = p.atom_sequences['A']
            chain_a_seq = str(chain_a_seq_o)
            a_occurrence = chain_a_seq.find(human_sequence[:-2]) # not all chains have the whole sequence
            if a_occurrence == -1 and pdb_id == '3H7P':
                a_occurrence = chain_a_seq.find(human_sequence[:62] + 'R' + human_sequence[63:]) # 3H7P has the R63K mutation for the isopeptide bond
            if chain_a_seq.startswith('MQ'):
                for x in range(len(human_sequence[:-2])):
                    if chain_a_seq[x] != human_sequence[x]:
                        print('MISMATCH AT POSITION {0} IN CHAIN A: Found {1}, expected {2}'.format(x, chain_a_seq[x], human_sequence[x]))
            assert(a_occurrence != -1)

            chain_b_seq_o = p.atom_sequences['B']
            chain_b_seq = str(chain_b_seq_o)
            b_occurrence = chain_b_seq.find(human_sequence[:-1]) # not all chains have the whole sequence
            if b_occurrence == -1 and pdb_id == '2XK5':
                b_occurrence = chain_b_seq.find(human_sequence[:38] + 'Q' + human_sequence[39:]) # 2XK5 has this extra mutation
            if chain_b_seq.startswith('MQ'):
                for x in range(len(human_sequence[:-2])):
                    if chain_b_seq[x] != human_sequence[x]:
                        print('MISMATCH AT POSITION {0} IN CHAIN B: Found {1}, expected {2}'.format(x, chain_b_seq[x], human_sequence[x]))
            assert(b_occurrence != -1)

            resfile.extend(['NATRO', 'start'])
            p19_1_key = chain_a_seq_o.order[a_occurrence + 19 - 1]
            p19_1_res = chain_a_seq_o.sequence[p19_1_key]
            e24_1_key = chain_a_seq_o.order[a_occurrence + 24 - 1]
            e24_1_res = chain_a_seq_o.sequence[e24_1_key]
            a28_1_key = chain_a_seq_o.order[a_occurrence + 28 - 1]
            a28_1_res = chain_a_seq_o.sequence[a28_1_key]
            assert(p19_1_res.ResidueAA == 'P' and e24_1_res.ResidueAA == 'E' and a28_1_res.ResidueAA == 'A')
            resfile.extend([
                '{0} A PIKAA S'.format(p19_1_res.ResidueID.strip()),
                '{0} A PIKAA D'.format(e24_1_res.ResidueID.strip()),
                '{0} A PIKAA S'.format(a28_1_res.ResidueID.strip()),
            ])

            p19_2_key = chain_b_seq_o.order[b_occurrence + 19 - 1]
            p19_2_res = chain_b_seq_o.sequence[p19_2_key]
            e24_2_key = chain_b_seq_o.order[b_occurrence + 24 - 1]
            e24_2_res = chain_b_seq_o.sequence[e24_2_key]
            a28_2_key = chain_b_seq_o.order[b_occurrence + 28 - 1]
            a28_2_res = chain_b_seq_o.sequence[a28_2_key]
            assert(p19_2_res.ResidueAA == 'P' and e24_2_res.ResidueAA == 'E' and a28_2_res.ResidueAA == 'A')
            resfile.extend([
                '{0} B PIKAA S'.format(p19_2_res.ResidueID.strip()),
                '{0} B PIKAA D'.format(e24_2_res.ResidueID.strip()),
                '{0} B PIKAA S'.format(a28_2_res.ResidueID.strip()),
            ])
            if pdb_id == '2XK5':
                p39_2_key = chain_b_seq_o.order[b_occurrence + 39 - 1]
                p39_2_res = chain_b_seq_o.sequence[p39_2_key]
                assert(p39_2_res.ResidueAA == 'Q')
                resfile.extend([
                    '{0} B PIKAA D'.format(p39_2_res.ResidueID.strip()),
                ])
        else:

            # Write new file
            pruned_lines = []
            for l in get_file_lines(filepath):
                if l[:6] == 'HETATM':
                    hetatm_types.add(l[17:20])
                if l[:6].strip() in chain_record_types:
                    assert(l[21] == 'A')
                if l[:6] != 'HETATM' or (l[17:20] == 'HOH'): # it seems fine to remove ZN ions in these cases
                    pruned_lines.append(l)
            write_file(pruned_filepath, '\n'.join(pruned_lines))
            p = PDB.from_filepath(pruned_filepath)
            assert(sorted(p.atom_sequences.keys()) == ['A'])

            assert(['A'] == p.atom_sequences.keys())
            chain_seq_o = p.atom_sequences['A']
            chain_sequence = str(p.atom_sequences['A'])
            first_occurrence = chain_sequence.find(human_sequence[1:])
            second_occurrence = chain_sequence[first_occurrence + 1:].find(human_sequence[:-2]) + (first_occurrence + 1)
            first_occurrence_seq = chain_sequence[first_occurrence:first_occurrence + len(human_sequence) - 1] # we subtract 1 because we threw the FME away
            linker_seq = chain_sequence[first_occurrence + len(human_sequence) - 1:second_occurrence]
            second_occurrence_seq = chain_sequence[second_occurrence:second_occurrence + len(human_sequence) - 2] # we subtract 2 because there are a couple of missing residues
            assert(first_occurrence_seq + linker_seq + second_occurrence_seq == chain_sequence)

            resfile.extend(['NATRO', 'start'])
            p19_1_key = chain_seq_o.order[first_occurrence + 19 - 2] # subtract an extra 1 for the missing MET
            p19_1_res = chain_seq_o.sequence[p19_1_key]
            e24_1_key = chain_seq_o.order[first_occurrence + 24 - 2]
            e24_1_res = chain_seq_o.sequence[e24_1_key]
            a28_1_key = chain_seq_o.order[first_occurrence + 28 - 2]
            a28_1_res = chain_seq_o.sequence[a28_1_key]
            assert(p19_1_res.ResidueAA == 'P' and e24_1_res.ResidueAA == 'E' and a28_1_res.ResidueAA == 'A')
            resfile.extend([
                '{0} A PIKAA S'.format(p19_1_res.ResidueID.strip()),
                '{0} A PIKAA D'.format(e24_1_res.ResidueID.strip()),
                '{0} A PIKAA S'.format(a28_1_res.ResidueID.strip()),
            ])

            p19_2_key = chain_seq_o.order[second_occurrence + 19 - 1]
            p19_2_res = chain_seq_o.sequence[p19_2_key]
            e24_2_key = chain_seq_o.order[second_occurrence + 24 - 1]
            e24_2_res = chain_seq_o.sequence[e24_2_key]
            a28_2_key = chain_seq_o.order[second_occurrence + 28 - 1]
            a28_2_res = chain_seq_o.sequence[a28_2_key]
            assert(p19_2_res.ResidueAA == 'P' and e24_2_res.ResidueAA == 'E' and a28_2_res.ResidueAA == 'A')
            resfile.extend([
                '{0} A PIKAA S'.format(p19_2_res.ResidueID.strip()),
                '{0} A PIKAA D'.format(e24_2_res.ResidueID.strip()),
                '{0} A PIKAA S'.format(a28_2_res.ResidueID.strip()),
            ])

        resfile_path = os.path.splitext(filepath)[0] + '.resfile'
        write_file(resfile_path, '\n'.join(resfile))
        print('\n'.join(resfile))

        score_filename = os.path.split(os.path.splitext(filepath)[0] + '.sc')[1]
        fixbb_script = '''
#!/bin/bash
#$ -cwd
#$ -r y
#$ -j n
#$ -l h_rt=24:00:00
#$ -l mem_free=2G
#$ -t 1
#$ -l arch=linux-x64

# This file was used to create the structures with the yeast sequence from the original human sequences.
# The command lines were provided by Samuel Thompson.

date
hostname

/netapp/home/shaneoconner/2015_q3_benchmarking/r58124/main/source/bin/fixbb.linuxgccrelease -s {0} -resfile {1} -ex1 -ex2 -extrachi_cutoff 0 -nstruct 1 -overwrite -linmem_ig 10 -minimize_sidechains -out:file:scorefile {2}

date
'''.format(os.path.split(pruned_filepath)[1], os.path.split(resfile_path)[1], score_filename)
        sub_script_path = os.path.splitext(filepath)[0] + '.sh'
        write_file(sub_script_path, fixbb_script)
        print('\n')


def add_headers():
    for pdb_id, details in year_2_cases.iteritems():
        source_filepath = details['filepath']
        target_filepath = details['filepath'].replace('.pdb', '.yeast.pdb')
        new_filepath = details['filepath'].replace('.pdb', '.yeast.headers.pdb')
        if read_file(target_filepath).strip().startswith('ATOM'):
            new_content = PDB.replace_headers(read_file(source_filepath), read_file(target_filepath))
            write_file(new_filepath, new_content)


def check_year_2_yeast_structures():
    d = copy.deepcopy(year_2_cases)
    for pdb_id, details in d.iteritems():
        filepath = details['filepath'].replace('.pdb', '.yeast.headers.pdb')
        details['db_filepath'] = filepath
        p = PDB.from_filepath(filepath)
        resfile = []
        if pdb_id != '2W9N':
            assert(p.atom_sequences.keys() == ['A', 'B'])
            a_sequence = str(p.atom_sequences['A'])
            if pdb_id == '3H7P':
                a_sequence = a_sequence[:62] + 'K' + a_sequence[63:]
            assert(a_sequence.find(yeast_sequence[:-2]) != -1)
            b_sequence = str(p.atom_sequences['B'])
            assert(b_sequence.find(yeast_sequence[:-2]) != -1)
        if pdb_id == '2W9N':
            assert(p.atom_sequences.keys() == ['A'])
            a_sequence = str(p.atom_sequences['A'])
            idx = a_sequence.find(yeast_sequence[1:])
            assert(idx != -1)
            idx = a_sequence[idx + 1:].find(yeast_sequence[:-2])
            assert(idx > 70)
            assert(len(a_sequence) < 160)
    return d


def import_year_2_structures():
    ppi_api = get_ppi_interface(read_file('pw'),
            #rosetta_scripts_path =  '/home/oconchus/t14benchmarking/r57934/main/source/bin/rosetta_scripts.linuxgccrelease',
            #rosetta_database_path = '/home/oconchus/t14benchmarking/r57934/main/database'
            )
    cases = check_year_2_yeast_structures()
    for pdb_id, details in sorted(cases.iteritems()):
        colortext.message('Importing {0}'.format(pdb_id))
        pprint.pprint(details)
        if pdb_id not in ['4S22', '3NOB']:
            continue
        # The UniProt mapping expects the reference PDB file to exist in the database so we will add it
        results = DDGdb.execute_select('SELECT * FROM PDBFile WHERE ID=%s', parameters=(pdb_id,))
        if len(results) == 0:
            fname = write_temp_file('/tmp', retrieve_pdb(pdb_id), suffix = '.pdb')
            print(fname)
            ppi_api.add_PDB_to_database(
                filepath = fname,
                pdbID = pdb_id,
                force = True,
                file_source = 'RCSB',
            )
            os.remove(fname)

        ppi_api.add_PDB_to_database(
            filepath = details['db_filepath'],
            pdbID = 'y' + pdb_id,
            derived_from = pdb_id,
            force = True,
            file_source = 'Rosetta',
            notes = "Created for the PUBS class at UCSF. Contact David Mavor, Kyle Barlow, Samuel Thompson, or Shane O'Connor for more details. This file is derived from an RCSB human ubiquitin structure but has been altered using the Rosetta fixbb application to use the yeast ubiquitin sequence. HETATM records and chains may also have been removed.",
            allow_missing_molecules = True
        )


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
        if pdb_id == '2W9N':
            continue

        existing_records = DDGdb.execute_select('''
            SELECT ID
            FROM PPIPDBSet
            INNER JOIN PPIPDBPartnerChain
            WHERE PPIPDBSet.PPComplexID = PPIPDBPartnerChain.PPComplexID
            AND PPIPDBSet.SetNumber = PPIPDBPartnerChain.SetNumber
            AND PPIPDBPartnerChain.PDBFileID=%s''', parameters=('y' + pdb_id,))
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
            PDBFileID = 'y' + pdb_id,
            Chain = 'A',
            NMRModel = None,
        )
        ppi_pdb_set_rpartner = dict(
            PPComplexID = complex_id,
            SetNumber = set_number,
            Side = 'R',
            ChainIndex = 0,
            PDBFileID = 'y' + pdb_id,
            Chain = 'B',
            NMRModel = None,
        )

        DDGdb.insertDictIfNew('PPIPDBSet', ppi_pdb_set, ['PPComplexID', 'SetNumber'])
        DDGdb.insertDictIfNew('PPIPDBPartnerChain', ppi_pdb_set_lpartner, ['PPComplexID', 'SetNumber', 'Side', 'ChainIndex'])
        DDGdb.insertDictIfNew('PPIPDBPartnerChain', ppi_pdb_set_rpartner, ['PPComplexID', 'SetNumber', 'Side', 'ChainIndex'])


def create_year_2_mutagenesis_records():
    ubiquitin_sequence = yeast_sequence

    for pdb_id, details in sorted(year_2_cases.iteritems()):
        if pdb_id == '2W9N':
            continue

        colortext.message(pdb_id)
        p = PDB(DDGdb.execute_select('SELECT Content FROM PDBFile WHERE ID=%s', parameters=('y' + pdb_id,))[0]['Content'])

        seq_a = str(p.atom_sequences['A'])
        seq_b = str(p.atom_sequences['B'])
        print('A:' + seq_a)
        print('B:' + seq_b)
        assert(seq_a.startswith('MQIF') or seq_a.startswith('SHMQIF'))
        assert(seq_a.endswith('VLRLRGG') or seq_a.endswith('VLRLRG') or seq_a.endswith('VLRLR') or seq_a.endswith('VLRL'))
        assert(seq_b.startswith('MQIF') or seq_b.startswith('SHMQIF'))
        assert(seq_b.endswith('VLRLRGG') or seq_b.endswith('VLRLRG') or seq_b.endswith('VLRLR') or seq_b.endswith('VLRL'))

        residue_ids_to_ignore = []
        typed_residue_ids_to_ignore = {'-1' : 'S', '0' : 'H'}
        for cr in details['conserved']:
            if pdb_id == '3H7P':
                residue_ids_to_ignore = [cr['ResidueID'] for cr in details['conserved']]
            else:
                typed_residue_ids_to_ignore[cr['ResidueID']] = 'K'

        mutageneses = p.generate_all_paired_mutations_for_position(['A', 'B'], residue_ids_to_ignore = residue_ids_to_ignore, typed_residue_ids_to_ignore = typed_residue_ids_to_ignore, silent = False)
        print('Number of mutations: {0}\n'.format(len(mutageneses)))

        # Get complex/PDB set IDs for the complex
        db_pdb_id = 'y' + pdb_id
        if pdb_id == '3H7P':
            complex_id = DDGdb.get_unique_record('SELECT ID FROM PPComplex WHERE LName=%s AND RName=%s', parameters=('Ubiquitin (yeast) K63R', 'Ubiquitin (yeast)'))['ID']
        else:
            complex_id = DDGdb.get_unique_record('SELECT ID FROM PPComplex WHERE LName=%s AND RName=%s', parameters=('Ubiquitin (yeast)', 'Ubiquitin (yeast)'))['ID']
        set_number = DDGdb.get_unique_record('SELECT DISTINCT SetNumber FROM PPIPDBPartnerChain WHERE PPComplexID=%s AND PDBFileID=%s AND Chain="A"', parameters=(complex_id, db_pdb_id))['SetNumber']

        for mutagenesis in mutageneses:
            id_string = ['y' + pdb_id]
            for mutation in sorted(mutagenesis):
                print(mutation.__dict__)
                print('{Chain}:{WildTypeAA}{ResidueID}{MutantAA}'.format(**mutation.__dict__))
                id_string.append('{Chain}:{WildTypeAA}{ResidueID}{MutantAA}'.format(**mutation.__dict__))
            id_string = '_'.join(id_string)

            # Use a transaction to prevent a partial deletion
            DDGdb._get_connection()
            con = DDGdb.connection
            try:
                with con:
                    cur = con.cursor()

                    # Add the mutagenesis record
                    pp_mutagenesis_record = dict(
                        PPComplexID = complex_id,
                        SKEMPI_KEY = id_string
                    )
                    pp_mutagenesis_id = DDGdb.transaction_insert_dict_auto_inc(cur, 'PPMutagenesis', pp_mutagenesis_record, unique_id_fields = ['SKEMPI_KEY'], check_existing = True)

                    #sql, params, record_exists = DDGdb.create_insert_dict_string('PPMutagenesis', pp_mutagenesis_record, PKfields = ['SKEMPI_KEY'], check_existing = True)
                    #if not record_exists:
                    #    cur.execute(sql, params)
                    #pp_mutagenesis_id = DDGdb.get_unique_record('SELECT ID FROM PPMutagenesis WHERE SKEMPI_KEY=%s', parameters=(id_string,))['ID']

                    pp_mutagenesis_mutations = []
                    for mutation in sorted(mutagenesis):
                        record_key = '{Chain} {WildTypeAA}{ResidueID}{MutantAA}'.format(**mutation.__dict__)
                        pp_mutagenesis_mutation = dict(
                            PPMutagenesisID = pp_mutagenesis_id,
                            RecordKey = '{Chain} {WildTypeAA}{ResidueID}{MutantAA}'.format(**mutation.__dict__),
                            ProteinID = None,
                            ResidueIndex = None,
                            WildTypeAA = mutation.WildTypeAA,
                            MutantAA = mutation.MutantAA,
                        )
                        pp_mutagenesis_mutation_id = DDGdb.transaction_insert_dict_auto_inc(cur, 'PPMutagenesisMutation', pp_mutagenesis_mutation, unique_id_fields = ['PPMutagenesisID', 'RecordKey'], check_existing = True)

                        pp_mutagenesis_pdb_mutation = dict(
                            PPMutagenesisID = pp_mutagenesis_id,
                            PPMutagenesisMutationID = pp_mutagenesis_mutation_id,
                            PPComplexID = complex_id,
                            SetNumber = set_number,
                            PDBFileID = db_pdb_id,
                            Chain = mutation.Chain,
                            WildTypeAA = mutation.WildTypeAA,
                            ResidueID = PDB.ResidueID2String(mutation.ResidueID),
                            MutantAA = mutation.MutantAA,
                        )
                        pp_mutagenesis_pdb_mutation_id = DDGdb.transaction_insert_dict_auto_inc(cur, 'PPMutagenesisPDBMutation', pp_mutagenesis_pdb_mutation, unique_id_fields = ['PPMutagenesisID', 'PDBFileID', 'Chain', 'ResidueID'], check_existing = True)

                        #sql, params, record_exists = DDGdb.create_insert_dict_string('PPMutagenesisMutation', pp_mutagenesis_mutation, PKfields = ['PPMutagenesisID', 'RecordKey'], check_existing = True)
                        #if not record_exists:
                        #    cur.execute(sql, params)
                        #pp_mutagenesis_mutation_id = DDGdb.get_unique_record('SELECT ID FROM PPMutagenesisMutation WHERE PPMutagenesisID=%s AND RecordKey=%s', parameters=(id_string,))['ID']

                        #print(
                        #id_string.append('{Chain}:{WildTypeAA}{ResidueID}{MutantAA}'.format(**mutation.__dict__))

            except Exception, e:
                raise colortext.Exception('An exception occurred removing the PredictionSet from the database: "%s".\n%s' % (str(e), traceback.format_exc()))


def check_counts():
    for pdb_id, details in sorted(year_2_cases.iteritems()):
        colortext.message(pdb_id)
        colortext.warning(len(DDGdb.execute_select('SELECT DISTINCT PPMutagenesisID FROM `PPMutagenesisPDBMutation` WHERE `PDBFileID`=%s', parameters = ('y' + pdb_id,))))


def create_year_2_user_dataset_records():
    user_dataset_name = 'DiPUBS'
    tnow = datetime.datetime.now()
    DDGdb.insertDictIfNew('UserDataSet', dict(
        TextID = user_dataset_name,
        UserID = 'oconchus',
        Description = 'A dataset for the ubiquitin runs for year 2 of the PUBS course at UCSF.',
        DatasetType = 'Binding affinity',
        FirstCreated = tnow,
        LastModified = tnow,

    ), ['TextID'])
    userdataset_id = DDGdb.get_unique_record('SELECT ID FROM UserDataSet WHERE TextID=%s', parameters=(user_dataset_name,))['ID']

    for pdb_id, details in sorted(year_2_cases.iteritems()):

        if pdb_id == '2W9N':
            continue

        # Get complex/PDB set IDs for the complex
        db_pdb_id = 'y' + pdb_id
        colortext.message(db_pdb_id)
        if pdb_id == '3H7P':
            complex_id = DDGdb.get_unique_record('SELECT ID FROM PPComplex WHERE LName=%s AND RName=%s', parameters=('Ubiquitin (yeast) K63R', 'Ubiquitin (yeast)'))['ID']
        else:
            complex_id = DDGdb.get_unique_record('SELECT ID FROM PPComplex WHERE LName=%s AND RName=%s', parameters=('Ubiquitin (yeast)', 'Ubiquitin (yeast)'))['ID']
        set_number = DDGdb.get_unique_record('SELECT DISTINCT SetNumber FROM PPIPDBPartnerChain WHERE PPComplexID=%s AND PDBFileID=%s AND Chain="A"', parameters=(complex_id, db_pdb_id))['SetNumber']

        for r in DDGdb.execute_select('SELECT DISTINCT PPMutagenesisID FROM `PPMutagenesisPDBMutation` WHERE `PDBFileID`=%s', parameters = (db_pdb_id,)):
            d = dict(
                UserDataSetID = userdataset_id,
                PPMutagenesisID = r['PPMutagenesisID'],
                PDBFileID = db_pdb_id,
                PPComplexID = complex_id,
                SetNumber = set_number,
                IsComplex = 1,
            )
            DDGdb.insertDictIfNew('UserPPDataSetExperiment', d, ['UserDataSetID', 'PPMutagenesisID', 'PDBFileID', 'PPComplexID', 'SetNumber'])
        print(len(DDGdb.execute_select('SELECT * FROM UserPPDataSetExperiment WHERE UserDataSetID=%s AND PDBFileID=%s', parameters=(userdataset_id, db_pdb_id))))


def add_year_2_prediction_set():

    ppi_api = get_ppi_interface(read_file('pw'),
        rosetta_scripts_path =  '/home/oconchus/2015_q3_benchmarking/r58124/main/source/bin/rosetta_scripts.linuxgccrelease',
        rosetta_database_path = '/home/oconchus/2015_q3_benchmarking/r58124/main/database'
        )

    # Add PredictionSet record
    prediction_set_id = 'DiPUBS: Complexes #1'
    ppi_api.add_prediction_set(prediction_set_id, description = 'DiPUBS: First run of the five complexes cases', halted = True, priority = 5, batch_size = 40, allow_existing_prediction_set = True)

    user_dataset_name = 'DiPUBS'
    ppi_api.add_prediction_run(prediction_set_id, user_dataset_name, extra_rosetta_command_flags = '-ignore_zero_occupancy false -ignore_unrecognized_res', protocol_id = None, tagged_subset = None, keep_hetatm_lines = False, test_run_first = False, show_full_errors = True)


def test_prediction_set():
    counts = {}
    ppi_api = get_ppi_interface(read_file('pw'))
    for j in ppi_api.get_queued_jobs('DiPUBS: Complexes #1', order_by = 'Cost', order_order_asc = False, include_files = True, truncate_content = None):
        counts[j['PDBFileID']] = counts.get(j['PDBFileID'], 0)
        counts[j['PDBFileID']] += 1
    pprint.pprint(counts)


test_prediction_set()

sys.exit(0)
if False:

    # I have (hopefully) written these functions so that they can be run multiple times without consequences e.g.
    #  - if you accidentally add the same mutation string a second time, it will be handled gracefully;
    #  - if you run add_predictions_by_pdb_id multiple times, it won't force jobs to be re-run.

    # Here is an example usage:

    ### Computational stage ###

    # Step 1: Add a new PDB
    ddG_connection.add_PDB_to_database('/kortemmelab/home/oconchus/ubiquitin/CUE.pdb', 'ub_CUE', file_source = 'Biosensor project', techniques = 'Rosetta model')
    ddG_connection.add_PDB_to_database('/kortemmelab/home/oconchus/ubiquitin/OTU.pdb', 'ub_OTU', file_source = 'Biosensor project', techniques = 'Rosetta model')
    ddG_connection.add_PDB_to_database('/kortemmelab/home/oconchus/ubiquitin/RPN13.pdb', 'ub_RPN13', file_source = 'Biosensor project', techniques = 'Rosetta model')
    ddG_connection.add_PDB_to_database('/kortemmelab/home/oconchus/ubiquitin/SH3.pdb', 'ub_SH3', file_source = 'Biosensor project', techniques = 'Rosetta model')
    ddG_connection.add_PDB_to_database('/kortemmelab/home/oconchus/ubiquitin/UQ_con.pdb', 'ub_UQcon', file_source = 'Biosensor project', techniques = 'Rosetta model')

    ddG_connection.add_PDB_to_database('/kortemmelab/home/oconchus/ubiquitin/1UBQ_yeast.pdb', 'uby_1UBQ', file_source = 'Biosensor project', techniques = 'Rosetta model')
    ddG_connection.add_PDB_to_database('/kortemmelab/home/oconchus/ubiquitin/OTU_yeast.pdb', 'uby_OTU', file_source = 'Biosensor project', techniques = 'Rosetta model')
    ddG_connection.add_PDB_to_database('/kortemmelab/home/oconchus/ubiquitin/RPN13_yeast.pdb', 'uby_RPN13', file_source = 'Biosensor project', techniques = 'Rosetta model')
    ddG_connection.add_PDB_to_database('/kortemmelab/home/oconchus/ubiquitin/SH3_yeast.pdb', 'uby_SH3', file_source = 'Biosensor project', techniques = 'Rosetta model')
    ddG_connection.add_PDB_to_database('/kortemmelab/home/oconchus/ubiquitin/UQ_con_yeast.pdb', 'uby_UQcon', file_source = 'Biosensor project', techniques = 'Rosetta model')
    ddG_connection.add_PDB_to_database('/kortemmelab/home/oconchus/ubiquitin/CUE_yeast.pdb', 'uby_CUE', file_source = 'Biosensor project', techniques = 'Rosetta model')

    # Step 2: Add a list of mutations for each PDB
    #
    # generate_all_point_mutations generates all possible point mutations for the specified chain.

    ubq_sequence = 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRG'
    for uc in ubiquitin_chains:
        # Generate a list of all possible mutations for the ubiquitin chains
        pdb_id = uc[0]
        chain_id = uc[1]
        p = PDB(ddG_connection.ddGDB.execute_select('SELECT Content FROM PDBFile WHERE ID=%s', parameters=(pdb_id,))[0]['Content'])
        #print('\n'.join(p.lines))
        assert(str(p.atom_sequences[chain_id]).find('MQIFV') != -1) and (str(p.atom_sequences[chain_id]).find('HLVLRLRG') != -1)
        mutant_list = p.generate_all_point_mutations_for_chain(chain_id)
        colortext.message('Generated %d mutations for %s.' % (len(mutant_list), pdb_id))
        colortext.warning('Adding mutants...')
        count = 0
        for mutant_mutation in mutant_list:
            ddG_connection.add_mutant(pdb_id, [mutant_mutation])
            if count % 20 == 0:
                sys.stdout.write('.')
            count += 1
        sys.stdout.write('\n')

if False:

    import pickle
    results = ddG_connection.ddGDB.execute_select('SELECT ID, PredictionSet, InputFiles FROM Prediction WHERE PredictionSet LIKE "Ubiquitin scan%%"', parameters=())
    for r in results:
        mutfile = pickle.loads(r['InputFiles'])['MUTFILE']
        if mutfile.find('total 1') == -1:
            print(mutfile)
            print('DELETE FROM Prediction WHERE ID=%s' %(r['ID']))
            #results = ddG_connection.ddGDB.execute('DELETE FROM Prediction WHERE ID=%s', parameters=(r['ID']))
            sys.exit(0)

if False:
    # Step 3:
    # This function will queue any mutations which have not previously run on the cluster or which are not in the queue.
    # The second parameter is a 'prediction set'. This is just a label used to group a set of predictions together. It makes sense
    # to group all predictions using the same protocol together.
    # The third parameter is the name of the DDG protocol. All protocols currently in the database are variants of protocol 16
    # from the Kellogg, Leaver-Fay, Baker paper (doi:10.1002/prot.22921).

    priority = 9 + len(ubiquitin_chains)
    count = 0
    for uc in ubiquitin_chains:
        priority -= 1 # assign different priorities to the different prediction sets
        ddG_connection.add_predictions_by_pdb_id(uc[0], uc[2], 'Protocol16 3.5.1 (talaris2013sc)', status = 'halted', priority = priority, KeepHETATMLines = False, strip_other_chains = False)
        count += 1
        #if count == 2:
        #    sys.exit(1)

if False:
    # Steps 2 and 3 can be repeated as often as you need however it is best to add as many mutations in step 2 first as this
    # will result in a better use of the cluster (larger array jobs).

    ### Analysis stage ###

    # To see the results quickly, you can use the get_flattened_prediction_results function e.g.
    colortext.message('Retrieving results')
    ddG_connection.get_flattened_prediction_results('FPP biosensor: protocol 16')

    # The create_abacus_graph_for_a_single_structure function will create a sorted graph of the results showing which sets of mutations are predicted to be more stable
    colortext.message('Analyzing results')
    ddG_connection.create_abacus_graph_for_a_single_structure('FPP biosensor: protocol 16', 'kellogg', 'total', graph_filename = 'L87Y_removed.png')

    # This is a simple test function which prints out the sequences of the monomer wildtype sequence and mutant sequence,
    # highlighting where they differ in case there was a mistake in the pipeline.
    # First, the output files are downloaded and extracted to the directory specified in the first argument.
    colortext.message('Testing results')
    ddG_connection.test_results('random_output_data', 'FPP biosensor: protocol 16')

    # This function creates a PyMOL session
    # The output files are downloaded and extracted to the directory specified in the second argument.
    # This function takes in a Prediction ID. These can be retrieved using the get_flattened_prediction_results function above.
    # The fourth argument is a task ID e.g. Protocol 16 generates 50 pairs of wildtype and mutant models by default, numbered
    # 1 to 50. This argument picks one of those pairs and creates a PyMOL session with aligned structures and with the mutations
    # highlighted.
    colortext.message('Creating PyMOL session')
    ddG_connection.create_pymol_session('test.pse', 'random_output_data', 53858, 25)

