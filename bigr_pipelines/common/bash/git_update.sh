# This script produces a bash script to create, update, and close issues

function create_gitlab_issue() {
    local URL="https://gitlab.com/bioinfo_gustaveroussy/bigr/${1}"
    local LABEL="${2}"
    local TITLE="${3}"
    local MESSAGE="${4}"

    CMD="glab issue open --title \"${TITLE}\" --description \"${DESCRIPTION}\" --yes --no-editor --label \"${LABEL}\""
    
    if [ ! -z ${GITLAB_USER:-} ]; then
        CMD+=" --assignee ${GITLAB_USER}"
    fi

    message CMD "${CMD}"
    # eval ${CMD}
}



function log_out() {
    local LOG_FILE="${1}"

    # TODO: Open issue
}

function log_error() {
    local LOG_FILE="${1}"

    # TODO: If no error, close issue with ok message
    # TODO: If error, close issue with error name
}


find logs/slurm/ -type f -name "*.log" | sort -r | while read LOG_FILE; do
    if [ "${LOG_FILE}" ~ *.err ]; then
        log_error "${LOG_FILE}"
    elif [ "${LOG_FILE}" ~ *.err ]; then
        log_out "${LOG_FILE}"
    else
        message WARNING "Dont' know how to handle ${FILE}"
    fi
done