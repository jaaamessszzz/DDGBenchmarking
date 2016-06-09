import os
import ast

def output_compiler():
    output_dir = '/kortemmelab/home/james.lucas/Structural_metrics_Outfiles'

    Assumed_set = set()
    Actual_set = set()

    for predID in range(94009, 95248):
        Assumed_set.add(predID)

    output_dict = {}
    fail_list = []

    for textfile in os.listdir(output_dir):
        if '.txt' in textfile:
            predID_txt = textfile[:-4]
            print predID_txt
            dict_entry = open(os.path.join(output_dir, textfile)).read()
            if 'Failure' in dict_entry:
                fail_list.append(int(ast.literal_eval(dict_entry + '}').keys()[0]))
                Actual_set.add(int(textfile[:-4]))
            else:
                print ast.literal_eval(dict_entry)
                output_dict[predID_txt] = ast.literal_eval(dict_entry)[predID_txt]
                Actual_set.add(int(textfile[:-4]))

    with open(os.path.join(output_dir, 'Structural_Metrics.py'), 'a') as final_out:
        final_out.write('RMSD_dict = %s' %str(output_dict))

    print '\nStuck Jobs'
    print Assumed_set - Actual_set

    print '\nCompleted Jobs (Even if they failed)'
    print Assumed_set & Actual_set

    print '\nFAIL.'
    print fail_list

output_compiler()