#!/usr/bin/env bash
set -euo pipefail

cd /data/vashah/seurat_explorer

export PATH=/opt/R/4.4.2/bin:$PATH
export R_LIBS_USER=/opt/R/library/4.4
export R_LIBS=/opt/R/library/4.4:/opt/R/4.4.2/lib64/R/library

unset RENV_PATHS_ROOT || true
unset RENV_PATHS_CACHE || true
unset RENV_CONFIG_SANDBOX_ENABLED || true
unset RENV_PROJECT || true

exec R --vanilla -e "
.libPaths(c('/opt/R/library/4.4', '/opt/R/4.4.2/lib64/R/library'));
print(.libPaths());
shiny::runApp('.', host='0.0.0.0', port=3838)
"
