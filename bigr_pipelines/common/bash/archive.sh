#!/usr/bin/env bash
set -e

# This adds several functions and variable in the environment
PIPELINE_PREFIX=$(readlink -e "$(dirname ${0})")
source "${PIPELINE_PREFIX}/messages.sh"
source "${PIPELINE_PREFIX}/environment.sh"

function move_to_archive() {
    local SOURCE="${1}"
    local DEST="${2}"

    message INFO "All data present in ${SOURCE} will be saved in ${DEST}"
    message INFO "This may take a while ..."
    COMMAND="rsync --exclude='.fq' --exclude='.fq.gz' --exclude='.fastq' --exclude='.fastq.gz' --exclude='logs' --exclude='.snakemake' --exclude='tmp' --exclude='log' --links --perms --times --group --omit-dir-times --verbose --checksum --recursive --update --progress --human-readable --partial ${SOURCE} ${DEST}"
    message CMD "${COMMAND}"
    eval ${COMMAND}
    #COMMAND="gh issue comment ${ISSUE_URL} --body \"${SOURCE} was saved at ${DEST}\""
    #message CMD "${COMMAND}"
    #eval ${COMMAND}
}

function remove_if_exists() {
    local TARGET="${1}"
    
    if [[ ! -z $(find -type d -name "${TARGET}") ]]; then
        message INFO "All data present in ${TARGET} directories will be deleted."
        COMMAND="find -type d -name ${TARGET} | while read SNAKEMAKE_TMP; do rm --recursive --force --verbose ${SNAKEMAKE_TMP:?}; done"
        message CMD "${COMMAND}"
        eval ${COMMAND}
    fi
}

PROJECT_PATH=$(readlink -e "${PWD}")
PROJECT_NAME=$(basename "${PROJECT_PATH}")

message WARNING "This script is interactive, do not leave."
# message INFO "I need your gitlab creditials to open, fill and close an issue about the project status."
# read -p "BiGR Gitlab identifier: " BIGR_GITLAB_ID
# read -s -p "BiGR Gitlab password: " BIGR_GITLAB_PW
# message INFO "We are working on the project ${PROJECT_NAME}"

DESTINATION="/mnt/isilon/data_bioinfo/Projets/bioinfo_core/${PROJECT_NAME}"
message INFO "By default, your project will be saved there: ${DESTINATION}"
message WARNING "If you do not belong to the Bioinformatics Core Facility, please change this destination."

while true; do
    read -p "New destination: (leave empty if you belong to BiGR) " newdest

    case "${newdest}" in
        "") message INFO "Destination remains default."; break;;
        * ) DESTINATION="${1}"; break;;
    esac
done

#COMMAND="gh issue create --assignee \"${BIGR_GITLAB_ID}\" --title \"Auto-archive\" --label \"In_Progress\" --body \"Archive of the project ${PROJECT_NAME} available at ${DESTINATION}\""
#message CMD "${COMMAND}"
#ISSUE_URL=$(eval ${COMMAND})

message INFO "Destination: ${DESTINATION}"
COMMAND="mkdir --parents --verbose \"${DESTINATION}\""
message CMD "${COMMAND}"
eval ${COMMAND}

move_to_archive "${PROJECT_PATH}/data_input/" "${DESTINATION}/data_input/"
move_to_archive "${PROJECT_PATH}/scripts/" "${DESTINATION}/scripts/"
move_to_archive "${PROJECT_PATH}/com/" "${DESTINATION}/com/"
move_to_archive "${PROJECT_PATH}/data_output" "${DESTINATION}/data_output/"
message INFO "You can now delete ${PROJECT_PATH}/data_input/ ${PROJECT_PATH}/script/ ${PROJECT_PATH}/com/ and ${PROJECT_PATH}/data_output safely."


while true; do
    read -p "Sould I remove directories named: .snakemake ? (y/n) " yn

    case "${yn}" in
        [yY] ) remove_if_exists ".snakemake"; break;;
        [nN] ) message INFO "Files have not been removed."; break;;
        * ) message ERROR "Unknown response.";;
    esac
done


while true; do
    read -p "Sould I remove directories named: tmp ? (y/n) " yn

    case "${yn}" in
        [yY] ) remove_if_exists "tmp"; break;;
        [nN] ) message INFO "Files have not been removed."; break;;
        * ) message ERROR "Unknown response.";;
    esac
done


while true; do
    read -p "Sould I remove directories named: logs ? (y/n) " yn

    case "${yn}" in
        [yY] ) remove_if_exists "logs"; break;;
        [nN] ) message INFO "Files have not been removed."; break;;
        * ) message ERROR "Unknown response.";;
    esac
done


while true; do
    read -p "Sould I remove directories named: log ? (y/n) " yn

    case "${yn}" in
        [yY] ) remove_if_exists "log"; break;;
        [nN] ) message INFO "Files have not been removed."; break;;
        * ) message ERROR "Unknown response.";;
    esac
done

#COMMAND="gh issue close ${ISSUE_URL} --comment 'Copy over.' --reason 'completed'"
# message CMD "${COMMAND}"
# eval ${COMMAND}

message INFO "Process over."