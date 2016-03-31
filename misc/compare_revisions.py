#!/usr/bin/python2

# The MIT License (MIT)
#
# Copyright (c) 2016 Shane O'Connor
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


import sys
import colortext


def read_file(filepath, binary = False):
    if binary:
        with open(filepath, 'rb') as f: contents = f.read()
    elif filepath.endswith('.gz'):
        f = gzip.open(filepath, 'r') # with...as fails for gzip.open in older versions of Python (e.g. v2.6.6)
        contents = f.read()
        f.close()
    else:
        with open(filepath, 'r') as f: contents = f.read()
    return contents


def parse_csv(filepath):
    '''pandas would be simpler to use here but just in case that is not installed on the system.'''
    lines = read_file(filepath).strip().split('\n')
    assert(lines)
    headers = lines[0].split(',')
    data = {}
    for l in lines[1:]:
        case = {}
        tokens = l.split(',')
        assert(len(tokens) == len(headers))
        for x in range(len(headers)):
            case[headers[x]] = tokens[x]
        case['DatasetID'] = int(case['DatasetID'])
        case['AbsoluteError_adj'] = float(case['AbsoluteError_adj'])
        case['Predicted_adj'] = float(case['Predicted_adj'])
        assert(case['AbsoluteError_adj'] >= 0)
        case['Experimental'] = float(case['Experimental'])
        data[case['DatasetID']] = case
    return data


def compare_mae(data1, data2, name1, name2):
    mae_diff = []
    case_constants = ['DatasetID', 'Experimental', 'HasGPMutation', 'MonomerLength', 'MutantAA', 'Mutations', 'NumberOfMutations', 'PDBFileID', 'PDBResolution', 'PDBResolutionBin', 'ResidueCharges', 'VolumeChange', 'WildTypeAA']
    assert(data1.keys() == data2.keys())
    for dataset_id, record1 in data1.iteritems():
        record2 = data2[dataset_id]
        for c in case_constants:
            if isinstance(record1[c], float):
                assert(abs(record1[c] - record2[c]) < 0.001)
            else:
                assert(record1[c] == record2[c])
        mae_diff.append( (record1['AbsoluteError_adj'] - record2['AbsoluteError_adj'], '{0} vs {1}'.format(str(round(record1['AbsoluteError_adj'], 2)), str(round(record2['AbsoluteError_adj'], 2))), dataset_id, record1['PDBFileID'], record1['Mutations']) )
    mae_diff = sorted(mae_diff)

    better_count = [0, 0]
    avg_count = [0.0, 0.0]
    print('{0} {1} {2} {3} {4} {5}'.format('Best method'.ljust(20), 'Delta'.ljust(6), 'Errors'.ljust(15), 'Case#'.ljust(6), 'PDB'.ljust(5), 'Mutations'))
    for m in mae_diff:
        best_method = '---'
        if m[0] < 0:
            best_method = name1
            better_count[0] += 1
            avg_count[0] += -m[0]
        elif m[0] > 0:
            best_method = name2
            better_count[1] += 1
            avg_count[1] += m[0]
        print('{0} {1} {2} {3} {4} {5}'.format(best_method[:20].ljust(20), str(round(m[0], 2)).ljust(6), str(m[1]).ljust(15), str(m[2]).ljust(6), m[3].ljust(5), m[4]))
    print('\n\n{0} performed better in {1} cases with average improvement {2}.'.format(name1, better_count[0], avg_count[0]/float(better_count[0])))
    print('{0} performed better in {1} cases with average improvement {2}.\n'.format(name2, better_count[1], avg_count[1]/float(better_count[1])))


if __name__ == '__main__':
    if not len(sys.argv) == 3 and not len(sys.argv) == 5:
        colortext.warning('\nThis script expects either:')
        colortext.warning('          - two filenames e.g. "{0} <file1.csv> <file2.csv>"'.format(sys.argv[0]))
        colortext.warning('          - or two filenames and two run names e.g. "{0} <file1.csv> <file2.csv> <run_name1> <run_name2>"'.format(sys.argv[0]))
        colortext.warning('e.g.')
        colortext.warning('python {0} july15v2_r58124_analysis_input.csv nov15v1_r58311_analysis_input.csv july15v2 nov15v1\n'.format(sys.argv[0]))
        sys.exit(1)

    name1 = 'Benchmark 1'
    name2 = 'Benchmark 2'
    if len(sys.argv) == 5:
        name1 = sys.argv[3]
        name2 = sys.argv[4]

    compare_mae(parse_csv(sys.argv[1]), parse_csv(sys.argv[2]), name1, name2)
