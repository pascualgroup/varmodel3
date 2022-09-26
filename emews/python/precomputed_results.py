from typing import Dict
import csv
import json

results: Dict[int, str] = {}


def load_results(results_file):
    with open(results_file) as fin:
        reader = csv.reader(fin)
        header = next(reader)
        instance_idx = header.index('instance')
        del header[instance_idx]
        for row in reader:
            instance = int(row[instance_idx])
            del row[instance_idx]
            result = {header[i]: float(v) for i, v in enumerate(row)}
            result['nbstrain'] = int(result['nbstrain'])
            result['nbgene'] = int(result['nbgene'])
            results[instance] = json.dumps(result)
    
    return results

def get_result(instance: int, results_file):
    global results
    if len(results) == 0:
        results = load_results(results_file)
    
    return results[instance]