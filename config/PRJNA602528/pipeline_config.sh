#!/usr/bin/env bash
: "${HF_HELIXFORGE_ROOT:?set HF_HELIXFORGE_ROOT}"
: "${HF_PACKAGE_ROOT:?set HF_PACKAGE_ROOT}"
export PROJECT_DIR="${HF_HELIXFORGE_ROOT}/pipelines/rnaseq"
export USER_SETTINGS_FILE="${HF_PACKAGE_ROOT}/config/PRJNA602528/user_settings.sh"
# shellcheck disable=SC1090
source "${HF_HELIXFORGE_ROOT}/pipelines/rnaseq/config/pipeline_config.sh"
