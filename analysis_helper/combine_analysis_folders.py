import sys
import os

output_dir = sys.argv[1]
assert( os.path.isdir(output_dir) )

subdirs = os.listdir( output_dir )

subdirs_data = []

if output_dir[-1] == '/':
    output_dir = output_dir[:-1]
run_name = os.path.split(output_dir)[-1]
    
for subdir in subdirs:
    s = subdir.split('_')
    n = None
    topx = None
    score_method = None
    for ss in s:
        if ss.startswith('n-'):
            n = int( ss[2:] )
        elif ss.startswith('topx-'):
            topx = int( ss[5:] )
        elif ss.startswith('method-'):
            score_method = int( ss[7:] )
            
    assert( n and topx and score_method )
    subdir_full_path = os.path.join(output_dir, subdir)
    metrics_file = os.path.join(subdir_full_path, '%s_metrics.txt' % run_name)
    assert( os.path.isfile(metrics_file) )
    pearsons_r = None
    mae = None
    fraction_correct = None
    with open(metrics_file, 'r') as f:
        for line in f:
            if line.startswith("Pearson's R"):
                pearsons_r = float( line.split(':')[1].strip().split()[0] )
            if line.startswith('MAE'):
                mae = float( line.split(':')[1].strip().split()[0] )
            if line.startswith('Fraction correct'):
                fraction_correct = float( line.split(':')[1].strip().split()[0] )
    assert( pearsons_r and mae and fraction_correct )
    subdirs_data.append( (mae, pearsons_r, fraction_correct, score_method, topx, subdir_full_path) )

print 'mae\tpearson_r\tfraction_correct\tscore_method\ttopx\tsubdir'
for mae, pearsons_r, fraction_correct, score_method, topx, subdir_full_path in sorted(subdirs_data):
    print '%.3f\t%.3f\t%.3f\t%d\t%d\t%s' % (mae, pearsons_r, fraction_correct, score_method, topx, subdir_full_path)


