#!/bin/bash
# run: ./00_process_run_in_hpc.sh $cfg_name $n_thread

#set -x
pid=$$
cfg_name=$1
if [[ $# -gt 1 ]]; then
	n_thread=$2
else
	n_thread=6;
fi
echo "Running MC-4C (python) pipeline for run [$cfg_name] using [$n_thread] threads for mapping."

jid_raw=$(qsub                    -terse -P compgen -N mc4c_${cfg_name}_01-cmb -l h_rt=40:00:00 -l h_vmem=10G -pe threaded 1 ~/bulk/bin/run_script.sh "python2 mc4c.py setReadIds ${cfg_name}")
echo "[1] Raw data processor is submitted with id: $jid_raw"

jid_spl=$(qsub -hold_jid $jid_raw -terse -P compgen -N mc4c_${cfg_name}_02-spl -l h_rt=40:00:00 -l h_vmem=10G -pe threaded 1 ~/bulk/bin/run_script.sh "python2 mc4c.py splitReads ${cfg_name}")
echo "[2] Read splitter is submitted with id: $jid_spl"

echo "Generating script for mapping fragments"
cmd_str=$(python2 mc4c.py mapFragments $cfg_name --return_command --n_thread $n_thread)
echo "submitting: $cmd_str"
jid_map=$(qsub -hold_jid $jid_spl -terse -P compgen -N mc4c_${cfg_name}_03-map -l h_rt=50:00:00 -l h_vmem=30G -pe threaded $n_thread ~/bulk/bin/run_script.sh "$cmd_str")
echo "[3] Fragment mapper is submitted with id: $jid_map"

jid_prc=$(qsub -hold_jid $jid_map -terse -P compgen -N mc4c_${cfg_name}_04-prc -l h_rt=50:00:00 -l h_vmem=30G -pe threaded 1 ~/bulk/bin/run_script.sh "python2 mc4c.py makeDataset ${cfg_name}")
echo "[4] Mapped fragment processor is submitted with id: $jid_prc"

jid_dup=$(qsub -hold_jid $jid_prc -terse -P compgen -N mc4c_${cfg_name}_05-dup -l h_rt=50:00:00 -l h_vmem=30G -pe threaded 1 ~/bulk/bin/run_script.sh "python2 mc4c.py removeDuplicates ${cfg_name}")
echo "[5] Duplicate filter is submitted with id: $jid_dup"

#### Plotting scripts
if [[ 1 -eq 1 ]]; then
echo "submitting jobs for plotting statistics"
    jid_stat01=$(qsub -hold_jid $jid_dup -terse -P compgen -N mc4c_${cfg_name}_07-01RdSizeDist -l h_rt=03:00:00 -l h_vmem=10G -pe threaded 1 ~/bulk/bin/run_script.sh "python2 mc4c.py getSumRep ${cfg_name}")
fi

#set +x