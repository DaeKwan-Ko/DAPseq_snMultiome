#!/usr/bin/env bash
set -euo pipefail

# 99_run_all_key_plotting.sh
# Run all key plotting scripts from the repository root.

Rscript scripts/08_key_plotting/01_plot_erf_dapseq_summaries.R
Rscript scripts/08_key_plotting/02_plot_evolutionary_summaries.R
Rscript scripts/08_key_plotting/03_plot_snmultiome_summaries.R
Rscript scripts/08_key_plotting/04_plot_network_summaries.R
Rscript scripts/08_key_plotting/05_prepare_source_data_tables.R
Rscript scripts/08_key_plotting/06_qc_check_plotting_inputs.R

echo "Key plotting workflow complete."
