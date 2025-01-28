import glob
import os

ROOT_DIR = "/mnt/isilon/sfgi/dudekm/"

MUNGE_LOGS = ROOT_DIR + "raw_data/GWAS/munged_sumstats"

def consolidate_results():
    annotations = [name for name in os.listdir("annotations") if os.path.isdir(f'annotations/{name}')]
    annotations.sort()
    print(annotations)

    # Format annotation name from directory name
    def get_annot_name(a):
        return a.replace('_new', '').replace('_overlap', '')

    outfile_traits = open('traits_summary.tsv', 'w')
    outfile_ldsr = open('ldsr_results.tsv', 'w')
    TRAITS_HEADER = ['Trait', 'Total_SNPs', 'Total_SNPs_after_munging', 'Regression_SNPs', 'Mean_Chi^2'] 
    LDSR_HEADER = ['Trait', 'Annotation', 'Prop_SNPs', 'Prop_h2', 'Enrichment', 'Enrichment_std_error', 'P', 'Coefficient', 'Coefficient_z']
    outfile_traits.write('\t'.join(TRAITS_HEADER) + '\n')
    outfile_ldsr.write('\t'.join(LDSR_HEADER) + '\n')

    munge_log_filenames = glob.glob(f'{MUNGE_LOGS}/*.log')
    munge_log_filenames.sort()

    for munge_log_filename in munge_log_filenames:
        trait = munge_log_filename.split('/')[-1].split('.')[0]
        munge_log_filename = f'{MUNGE_LOGS}/{trait}.log'

        print(f'Getting results for trait "{trait}"...', flush=True)

        # Read munge summary
        munge_log = open(munge_log_filename)
        for i, line in enumerate(munge_log):
            if 'SNPs from --sumstats file' in line:
                total_snps = line.split('Read ')[-1].split(' SNPs from --sumstats file')[0]
            elif 'Writing summary statistics' in line:
                total_snps_after_munging = line.split('SNPs (')[-1].split(' with nonmissing beta)')[0]
            elif 'Mean chi^2 = ' in line:
                mean_chi2 = line.strip().split('Mean chi^2 = ')[-1]
        munge_log.close()

        # Read ldsr results for each annotation
        regression_snps = 'NA'
        for annot in annotations:
            result_filename = f'annotations/{annot}/out/{trait}.{annot}.results'
            ldsr_log_filename = f'annotations/{annot}/out/{trait}.{annot}.log'
            if not os.path.isfile(result_filename):
                print(f'WARNING: could not find file {result_filename}', flush=True)
                p = 'NA'
                enrichment = 'NA'
            else:
                result = open(result_filename)
                ldsr_log = open(ldsr_log_filename)

                # We only need to get regression_snps once, so just save the last one
                for i, line in enumerate(ldsr_log):
                    if ' SNPs remain)' in line:
                        regression_snps = line.split('(')[-1].split(' SNPs remain)')[0]
                
                for i, line in enumerate(result):
                    if i == 1:
                        _, prop_snps, prop_h2, _, enrichment, enrichment_std_error, p, _, coefficient, coefficient_z = line.strip().split('\t')
                        break
                
            ldsr_row = [trait, get_annot_name(annot), prop_snps, prop_h2, enrichment, enrichment_std_error, p, coefficient, coefficient_z]
            outfile_ldsr.write('\t'.join(ldsr_row) + '\n')

            result.close()
            ldsr_log.close()

        trait_row = [trait, total_snps, total_snps_after_munging, regression_snps, mean_chi2]
        outfile_traits.write('\t'.join(trait_row) + '\n')
        

    outfile_traits.close()
    outfile_ldsr.close()


if __name__ == '__main__':
    consolidate_results()
