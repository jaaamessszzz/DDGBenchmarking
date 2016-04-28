import os
import prody
import re
import requests
import klab.bio.rcsb as rcsb
import klab.bio.alignment as align
import klab.bio.pdb as PDB
import pandas as pd
import Bio.PDB

#From Kyle's Finalize.py
def read_mutations_resfile(filenum_dir):
    resfile = os.path.join(filenum_dir, 'mutations_repack.resfile')
    mutations = []
    with open(resfile, 'r') as f:
        post_start = False
        for line in f:
            if post_start:
                line = line.strip()
                pdb_resnum, chain, pikaa, mut_res = line.split()
                mutations.append( [pdb_resnum, chain, pikaa, mut_res] )
            elif line.startswith('start'):
                post_start = True
    return mutations

def generate_mut_seq(pdbdir, hv):
    blast_dict = {}
    mutation_dict = {}
        
    mutations = read_mutations_resfile(pdbdir)
        
    for mutation in mutations:
        mutation_dict[mutation[0] + mutation[1]] = mutation[3]
    
    #Dictionary: 3- to 1-letter code
    res_dict = {
        'ALA':'A',
        'CYS':'C',
        'ASP':'D',
        'GLU':'E',
        'PHE':'F',
        'GLY':'G',
        'HIS':'H',
        'ILE':'I',
        'LYS':'K',
        'LEU':'L',
        'MET':'M',
        'ASN':'N',
        'PRO':'P',
        'GLN':'Q',
        'ARG':'R',
        'SER':'S',
        'THR':'T',
        'VAL':'V',
        'TRP':'W',
        'TYR':'Y'
    }
    
    for chain in hv:
        chain_string = ''
        for residue in chain:
            current_res =  str(residue.getResnum()) + str(residue.getChid())
            if current_res not in mutation_dict.keys():
                chain_string = chain_string + res_dict[residue.getResname()]
            else:
                chain_string = chain_string + mutation_dict[str(residue.getResnum()) + str(residue.getChid())]
        blast_dict[chain.getChid()] = chain_string
    
    return blast_dict

def get_PDBs(pdbdir):
    skip_blast = False
    
    for pdb_file in os.listdir(pdbdir):
        if pdb_file.endswith('.pdb'):
            print "\n***Now BLASTing %s***" % pdb_file
            actual_pdb_file = pdb_file
            output_dir = os.path.join('/kortemmelab/home/james.lucas/Mutant_PDBs', pdb_file[:-4])
            try:
                os.makedirs(output_dir)        
            except:
                print 'File directory %s already exists!!!' % pdb_file[:-4]
                skip_blast = True
            
            if skip_blast == False:
                pdb_chains = re.split(r'_', pdb_file)[1][:-4]
                pdb_id = re.split(r'_', pdb_file)[0]
                hv = prody.parsePDB(os.path.join(pdbdir, pdb_file)).getHierView()

                blast_dict = generate_mut_seq(pdbdir, hv)
                hit_list = []
                for sequence in blast_dict:
                    set_dict = {}
                    temp_set = set()
                    blast_me = prody.blastPDB(blast_dict[sequence], timeout = 600)
                    hits_dict = blast_me.getHits(percent_identity=95)
                    for hit in hits_dict:
                        temp_set.add(hit)

                    hit_list.append(temp_set)

                common_set = hit_list[0]

                for entry in hit_list[1:]:
                    common_set.intersection_update(entry)

                for hit in common_set:
                    rcsb.download_pdb(hit, output_dir)

                alignments(actual_pdb_file, common_set, pdbdir)
            else:
                print 'Blast skipped!'

def alignments(pdb_file, common_set, pdbdir):
    for pdb_id in common_set:
        clustal(pdb_file, pdb_id, pdbdir)

#Shane's ClustalO code
def clustal(pdb_WT, pdb_target, pdbdir):
    import pprint
    from klab.bio.clustalo import SequenceAligner
    from klab.bio.uniprot import UniParcEntry
    from klab.bio.basics import SubstitutionScore, PDBResidue
    from klab import colortext
    
    symbols = SubstitutionScore.clustal_symbols
    
    wt_pdb_path = os.path.join('/kortemmelab', 'home', 'james.lucas', 'jobs', '160322-james-backrub-rscript-full', 'data', pdbdir, 
                                            '{0}.pdb'.format(pdb_WT[:-4]))

    p1 = PDB.PDB.from_filepath(wt_pdb_path)
    p2 = PDB.PDB.from_filepath(os.path.join('/kortemmelab', 'home', 'james.lucas', 'Mutant_PDBs', pdb_WT[:-4], '{0}.pdb'.format(pdb_target)))

    sa = SequenceAligner()
    for chain_id, seq in (p1.seqres_sequences or p1.atom_sequences).iteritems():
        sa.add_sequence('{0}_{1}'.format(pdb_WT[:-4], chain_id), str(seq))
    for chain_id, seq in p2.seqres_sequences.iteritems():
        sa.add_sequence('{0}_{1}'.format(pdb_target, chain_id), str(seq))
    sa.align()
    
    # Prints the raw ClustalO output (for debugging)
    #colortext.message('\nClustalO output')
    #print(sa.alignment_output)

    # Find all differing residues
                    
    for chain_id_WT, seq_WT in (p1.seqres_sequences or p1.atom_sequences).iteritems():
        best_match_dict = sa.get_best_matches_by_id( '%s_%s' %(pdb_WT[:-4], chain_id_WT), cut_off = 90.0)
        
        for best_match_pdb in best_match_dict.keys():
            sa2 = SequenceAligner()
            sa2.add_sequence('{0}_{1}'.format(pdb_WT[:-4], chain_id_WT), str(seq_WT))
            sa2.add_sequence('{0}_{1}'.format(best_match_pdb[:-2], best_match_pdb[-1:]), str(p2.seqres_sequences[best_match_pdb[-1:]]))
            asddfasdf=p2.seqres_sequences[best_match_pdb[-1:]]
            sa2.align()
            resmap, clustal_matches = sa2.get_residue_mapping()
            
            #print clustal_matches
            raw_mut_list = []
            for first_idx, v in sorted(clustal_matches.iteritems()):
                # first_idx is an index into the first sequence, v is a SubstitutionScore object
                if symbols[v.clustal] != '*':
                    # get the index into the second sequence
                    second_idx = resmap[first_idx]

                    pdb_res1 = seq_WT.sequence[seq_WT.order[first_idx - 1]] # low-level access - there is probably a utility function to do this but I was hurrying
                    pdb_res2 = asddfasdf.sequence[asddfasdf.order[second_idx - 1]]
                    #print('  {0} -> {1}'.format(super(PDBResidue, pdb_res1).__repr__(), pdb_res2))
                    
                    raw_mut_list.append('  {0} -> {1}'.format(super(PDBResidue, pdb_res1).__repr__(), pdb_res2))
            mut_list_interpreter(raw_mut_list, pdb_WT, pdb_target)

#Takes list of mutation strings and converts it into mutations format
#If there is a 100% match in df, add remark in a new column with matched PDBID
def mut_list_interpreter(formatted_mut_list, pdb_WT, pdb_target):
    mut_string_temp = ''
    for element in formatted_mut_list:
        split_input = re.split(':| ', element.strip())
        residue_position = '    ' + split_input[1]
        mut_string_temp = mut_string_temp + '%s %s%s %s; ' %(split_input[0], split_input[2], residue_position[-4:], split_input[6])
    mut_string = mut_string_temp[:-2]
    
    output = open('/kortemmelab/home/james.lucas/Mutant_PDBs/Resfile_Mutations_Completed.txt', 'a')
    df = pd.read_csv('/kortemmelab/home/james.lucas/Mutant_PDBs/Resfile_Mutations.csv')
    for index, row in df.iterrows():
        if row[1].strip()== mut_string: 
            if pdb_WT[:4] == row[2].strip():
                output.write('%s,%s,%s,%s\n' % (row[0], pdb_WT, row[1], pdb_target) )
                print '%s,%s,%s,%s\n' % (row[0], pdb_WT, row[1], pdb_target)
            
def main():
    os.chdir('/kortemmelab/home/james.lucas/jobs/160322-james-backrub-rscript-full/data/')
    cwd = os.getcwd()
    output = open('/kortemmelab/home/james.lucas/Mutant_PDBs/Resfile_Mutations_Completed.txt', 'a')
    output.write('DatasetID,PDBFileID,Mutations,Matches\n')
    for pdbdir in os.listdir(cwd):
        if os.path.isdir(pdbdir):
            get_PDBs(pdbdir)
            
main()
