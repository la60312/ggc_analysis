import os
import argparse


parser = argparse.ArgumentParser(description='Run PEDIA analysis')
parser.add_argument('Train_path', help='path of training data')
parser.add_argument('Test_path', help='path of testing data')
parser.add_argument('-o', '--output', default=".", help='Path of output')
args = parser.parse_args()

train = args.Train_path
#train = 'pedia_train_jsons'
test_dir = args.Test_path
#test_dir = 'pedia_ggc_test_jsons'

result_dir = args.output

features = ['0_2', '3_4', '1_3_4', '1_2_3_4', '0_2_3_4', '0_1_3_4']
out_dirs = ['ggc_c_p_b/', 'ggc_f_c_g/','ggc_f_g/','ggc_f/','ggc_c/','ggc_g/']

for idx, f in enumerate(features):
    out_dir = os.path.join(result_dir, out_dirs[idx])
    cmd = 'python3 classifier/pedia.py {} 1KG -t {} -e {} -o {} -p 5'.format(train, test_dir, f, out_dir)
    print(cmd)
    os.system(cmd)

out_dir = os.path.join(result_dir, 'ggc_pedia')
cmd = 'python3 classifier/pedia.py {} 1KG -t {} -o {} -p 5'.format(train, test_dir,  out_dir)
print(cmd)
os.system(cmd)
