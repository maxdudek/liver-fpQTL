import glob
import os

ROOT_DIR = "/mnt/isilon/sfgi/dudekm/"

def consolidate_results(dir='.'):
    print(f'dir = {dir}')
    annotations = [name for name in os.listdir(f'{dir}/annotations') if os.path.isdir(f'{dir}/annotations/{name}')]
    annotations.sort()
    print(annotations)

    # Format annotation name from directory name
    def get_annot_name(a):
        return a

    outfile_ldsr = open(f'{dir}/ldsr_continuous_results.tsv', 'w')
    LDSR_HEADER = ['Trait', 'Annotation', 'Regression_Model', 'Coefficient', 'Coefficient_StdErr', 'Coefficient_z', 'Q5_prop_h2g', 'Q5_prop_h2g_se', 'Q5_enr', 'Q5_enr_se', 'Q5_enr_pval']
    outfile_ldsr.write('\t'.join(LDSR_HEADER) + '\n')

    for annot in annotations:
        print(f'\tAnnot = {annot}')
        
        for regression_model_dir in ['out_baseline', 'out_ocr', 'out_ocr_baseline']:
            regression_model = regression_model_dir.replace('out_', '')
            ldsr_out_files = glob.glob(f'{dir}/annotations/{annot}/{regression_model_dir}/*.results')
            for ldsr_out_file in ldsr_out_files:
                trait = ldsr_out_file.split('/')[-1].split(f'.{annot}.results')[0]
                q5_outfile_name = f'{dir}/annotations/{annot}/{regression_model_dir}/{trait}.{annot}.q5.txt'

                # print(f'\t\tGetting results for trait "{trait}"...', flush=True)

                result = open(ldsr_out_file)
                for i, line in enumerate(result):
                    if i == 1:
                        _, _, _, _, _, _, _, coefficient, coefficient_std_err, coefficient_z = line.strip().split('\t')
                        break
                result.close()

                N_QUANTILES = 5
                prop_h2g, prop_h2g_se, enr, enr_se, enr_pval = ([None]*N_QUANTILES,)*5
                q5_file = open(q5_outfile_name)
                for q, line in enumerate(q5_file):
                    if q > 0:
                        _, _, prop_h2g[q-1], prop_h2g_se[q-1], enr[q-1], enr_se[q-1], enr_pval[q-1] = line.strip().split('\t')
                    
                ldsr_row = [trait, get_annot_name(annot), regression_model, coefficient, coefficient_std_err, coefficient_z,
                            ';'.join(prop_h2g), ';'.join(prop_h2g_se), ';'.join(enr), ';'.join(enr_se), ';'.join(enr_pval)]
                outfile_ldsr.write('\t'.join(ldsr_row) + '\n')

                
    
    outfile_ldsr.close()

if __name__ == '__main__':
    consolidate_results()
