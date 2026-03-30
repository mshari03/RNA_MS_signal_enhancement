import os
import glob
import subprocess
from datetime import datetime

# all outputs are placed into a folder labeled with the date and time
current_datetime = datetime.now()
datetimestr = current_datetime.strftime("%Y-%m-%d-%H-%M-%S")
subprocess.run(f"mkdir outputs/{datetimestr}", shell=True)
subprocess.run(f"mkdir outputs/{datetimestr}/in_silico_output", shell=True)
subprocess.run(f"mkdir outputs/{datetimestr}/final_reports", shell=True)
subprocess.run(f"mkdir outputs/{datetimestr}/FDR_reports", shell=True)

#### PARAMETERS ####
# enter your data folder here
input_folder = "20250717"
theoretical_digest_file = "Digest_31mer"

# ---- in silico digestion (do this once per modified sequence) ----
enzyme_CL = f"pythonw in_silico_digestion/1_enzyme.py \
            --RNA_sequences inputs/{input_folder}/31mer.fasta \
            --enzyme T1 \
            --miss 3 \
            --nonspecific_min_length 3 \
            --nonspecific_max_length 12 \
            --cleaved_fragments_5end_chem OH \
            --cleaved_fragments_3end_chem OH P cP \
            --RNA_5end_chem OH \
            --RNA_3end_chem OH"
print(enzyme_CL)
subprocess.run(enzyme_CL, shell=True)
subprocess.run(f"cp output.1 outputs/{datetimestr}/in_silico_output/enzyme_output", shell=True)

modify_CL = f"pythonw in_silico_digestion/2_modify.py \
            --nts_light in_silico_digestion/nts_light.xlsx"
subprocess.run(modify_CL, shell=True)
subprocess.run(f"cp output.2 outputs/{datetimestr}/in_silico_output/modify_output", shell=True)

consolidate_CL = "pythonw in_silico_digestion/3_consolidate.py \
            --MS MS2\
            --min_length 3"
subprocess.run(consolidate_CL, shell=True)
subprocess.run(f"cp output.3.MS2 outputs/{datetimestr}/in_silico_output/consolidate_output", shell=True)
subprocess.run("pythonw in_silico_digestion/3.5_decoy.py", shell=True)
subprocess.run(f"cp output.3.MS2 outputs/{datetimestr}/in_silico_output/decoy_output", shell=True)

calc_mass_CL = "pythonw in_silico_digestion/4_calc_mass.py \
        --ion_mode - \
        --nts_light in_silico_digestion/nts_light.xlsx \
        --MS1_charges in_silico_digestion/charges_MS1.txt \
        --MS2_charges in_silico_digestion/charges_MS2.txt \
        --MS1_mzlow 400 \
        --MS1_mzhigh 3200 \
        --MS2_mzlow 100 \
        --MS2_mzhigh 3200 "
subprocess.run(calc_mass_CL, shell=True)

# file organization and deletion
subprocess.run(f"cp {theoretical_digest_file}.txt outputs/{datetimestr}/in_silico_output/{theoretical_digest_file}.txt", shell=True)
subprocess.run("rm output.1", shell=True)
subprocess.run("rm output.2", shell=True)
subprocess.run("rm output.3.MS2", shell=True)
subprocess.run(f"rm {theoretical_digest_file}.txt", shell=True)


# ---- matching scoring ----
# loop through this for every input file
for filepath in glob.glob(f"inputs/{input_folder}/*.mgf"):
        filename = os.path.splitext(os.path.basename(filepath))[0]
        matching_CL = f"pythonw matching_scoring/pytheas_matching.py \
        --theoretical_digest outputs/{datetimestr}/in_silico_output/{theoretical_digest_file}.txt \
        --mgf_file {filepath} \
        --isotopic_species light \
        --MS1_ppm 50 \
        --MS2_ppm 50 \
        --MS1_mz_minimum 400 \
        --MS1_mz_maximum 3200 \
        --MS2_mz_minimum 100 \
        --MS2_mz_maximum 3200 \
        --MS2_peak_int_min all \
        --MS2_peak_num_max all \
        --alpha 0 \
        --beta_increment 0.075 \
        --precursor_isotopologues y \
        --use_charges_mgf y \
        --FDR_isotopic_species light \
        --only_targets_with_decoys y \
        --sequence_lengths_FDR all"
        subprocess.run(matching_CL, shell=True)

        # file organization for matching and scoring
        subprocess.run(f"mkdir outputs/{datetimestr}/match_outputs_{filename}", shell=True)
        subprocess.run(f"cp match_output_{filename}.txt outputs/{datetimestr}/match_outputs_{filename}/match_output_{filename}.txt", shell=True)
        subprocess.run(f"cp decoys_{filename}.csv outputs/{datetimestr}/match_outputs_{filename}/decoys_{filename}.csv", shell=True)
        subprocess.run(f"cp targets_{filename}.csv outputs/{datetimestr}/match_outputs_{filename}/targets_{filename}.csv", shell=True)
        subprocess.run(f"cp log.txt outputs/{datetimestr}/match_outputs_{filename}/log_{filename}.txt", shell=True)
        subprocess.run(f"rm match_output_{filename}.txt", shell=True)
        subprocess.run(f"rm decoys_{filename}.csv", shell=True)
        subprocess.run(f"rm targets_{filename}.csv", shell=True)
        subprocess.run(f"rm log.txt", shell=True)

        # ---- final report ----
        report_CL = f"pythonw final_report/pytheas_final_report.py \
        --match_file outputs/{datetimestr}/match_outputs_{filename}/match_output_{filename}.txt \
        --visualize_decoys y \
        --modified_only n \
        --only_unique_positions n \
        --Sp_minimum_cutoff 0 \
        --dSp_maximum_cutoff 0 \
        --rank_max 99 \
        --remove_redundant_sequences_X n \
        --dSp2_minimum_cutoff 0 "
        subprocess.run(report_CL, shell=True)

        # file organization for final report
        subprocess.run(f"cp final_report_{filename}.csv outputs/{datetimestr}/final_reports/{filename}_final_report.csv", shell=True)
        subprocess.run(f"cp FDR_{filename}.csv outputs/{datetimestr}/FDR_reports/{filename}_FDR.csv", shell=True)
        subprocess.run(f"rm FDR_{filename}.csv", shell=True)
        subprocess.run(f"rm final_report_{filename}.csv", shell=True)
        
subprocess.run(f"rm seq_output", shell=True)

# For some reason these other pytheas scripts are throwing errors, need to debug them further for complete automation
# visualization
# visualization_CL = f"pythonw visualization_spectra/pytheas_visualization_html.py \
#     --digest_file outputs/{datetimestr}/in_silico_output/{theoretical_digest_file}.txt\
#     --mgf_file inputs/{input_folder}/{filename}.txt \
#     --match_file outputs/{datetimestr}/match_outputs_{filename}/match_output_{filename}.txt"
# subprocess.run(visualization_CL, shell=True)


# # sequence mapping
# mapping_CL = f"pythonw sequence_mapping/pytheas_mapping.py \
#     --nts_alphabet nts_light.csv \
#     --input_file outputs/{datetimestr}/{filename}_final_report.txt \
#     --input_sequences seq_output \
#     --min_length 5 \
#     --Sp_cutoff 0"
# subprocess.run(mapping_CL, shell=True)


# # other optional pytheas scripts could be added




