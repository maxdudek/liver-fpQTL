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

    outfile_ldsr = open(f'{dir}/ldsr_results.tsv', 'w')
    LDSR_HEADER = ['Trait', 'Annotation', 'Regression_Model', 'Prop_SNPs', 'Prop_h2', 'Enrichment', 'Enrichment_std_error', 'P', 'Coefficient', 'Coefficient_z']
    outfile_ldsr.write('\t'.join(LDSR_HEADER) + '\n')

    for annot in annotations:
        print(f'\tAnnot = {annot}')
        
        for regression_model_dir in ['out_baseline', 'out_ocr', 'out_ocr_baseline']:
            regression_model = regression_model_dir.replace('out_', '')
            ldsr_out_files = glob.glob(f'{dir}/annotations/{annot}/{regression_model_dir}/*.results')
            for ldsr_out_file in ldsr_out_files:
                trait = ldsr_out_file.split('/')[-1].split(f'.{annot}.results')[0]

                # print(f'\t\tGetting results for trait "{trait}"...', flush=True)

                result = open(ldsr_out_file)
                for i, line in enumerate(result):
                    if i == 1:
                        _, prop_snps, prop_h2, _, enrichment, enrichment_std_error, p, _, coefficient, coefficient_z = line.strip().split('\t')
                        break
                    
                ldsr_row = [trait, get_annot_name(annot), regression_model, prop_snps, prop_h2, enrichment, enrichment_std_error, p, coefficient, coefficient_z]
                outfile_ldsr.write('\t'.join(ldsr_row) + '\n')

                result.close()
    
    outfile_ldsr.close()

if __name__ == '__main__':
    consolidate_results()
