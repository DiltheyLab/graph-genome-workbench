## Parsing json
import json 
import sys
import argparse

def extract_data(data_dict):
    tp_base = data_dict['TP-base']
    tp_call = data_dict['TP-comp']
    fp = data_dict['FP']
    fn = data_dict['FN']
    precision = data_dict['precision']
    recall = data_dict['recall']
    fscore = data_dict['f1']
    return [None, tp_base, tp_call, fp, fn, precision, recall, fscore]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='parse_json_to_txt.py', description=__doc__)
    parser.add_argument('json_file', metavar='JSON', help='Name of JSON file')
    args = parser.parse_args()
    with open(args.json_file, 'r') as file:
        data = json.load(file)

    print('\t'.join(["Threshold", 'True-pos-baseline', 'True-pos-call', 'False-pos', 'False-neg', 'Precision', 'Sensitivity', 'F-measure']))
    print('----------------------------------------------------------------------------------------------------')
    print('\t'.join([str(v) for v in extract_data(data)]))


