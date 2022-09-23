#!/usr/bin/env bash
set -e

# This adds several functions and variable in the environment
PIPELINE_PREFIX=$(readlink -e "$(dirname ${0})")
source "${PIPELINE_PREFIX}/messages.sh"
source "${PIPELINE_PREFIX}/environment.sh"

MAX_THREADS=5

message WARNING "This script is interactive, do not leave."
message INFO "I need your password to access datasets..."
COMMAND="iinit"

while true; do
    read -p "Sould I run ${COMMAND} command ? (y/n) " yn

    case "${yn}" in
        [yY] ) message CMD "${COMMAND}"; eval ${COMMAND}; break;;
        [nN] ) message INFO "Skipping iinit."; break;;
        * ) message ERROR "Unknown response.";;
    esac
done



message INFO "I will use ${MAX_THREADS} threads to search information about the project ${1}"

COMMAND="imeta qu -C 'projectName' like \"${1}\" | grep -v \"\-\-\-\" | cut -f2 -d\" \" > datasetList"
message CMD "${COMMAND}"
eval ${COMMAND}
message INFO "I have found $(wc -l datasetList) datasets (including re-sequencing)"

while true; do
    read -p "Do you feel like its ok to continue ? (y/n) " yn

    case "${yn}" in
        [yY] ) message INFO "Ok, I shall proceeed"; break;;
        [nN] ) message ERROR "Exiting..."; exit ;;
        * ) message ERROR "Unknown response.";;
    esac
done


my_iquest() { iquest "%s/%s" "SELECT COLL_NAME,DATA_NAME WHERE COLL_NAME like '${1}%'"; }
export -f my_iquest
COMMAND="bash ${PIPELINE_PREFIX}/multiprocess_bash_command.sh -p ${MAX_THREADS} -f my_iquest $(cat datasetList | tr $'\n' ' ') > dataset_paths.txt"
message CMD "${COMMAND:0:50}..."
eval ${COMMAND}


while true; do
    read -p "Is this dataset paried ? (y/n) " yn

    case "${yn}" in
        [yY] ) message INFO "Ok, I will pair fastq files alphabetically"; break;;
        [nN] ) message ERROR "Nothing more to do."; exit ;;
        * ) message ERROR "Unknown response.";;
    esac
done


COMMAND='grep -vP ".*I[1,2]_001\..*\.gz" dataset_paths.txt | sort | uniq | paste - - > paired_dataset_paths.txt'
message CMD "${COMMAND}"
eval ${COMMAND}

while true; do
    read -p "Should I try to guess resequenced samples ? (y/n) " yn

    case "${yn}" in
        [yY] ) message INFO "Ok, I will guess multi-sequenced samples"; break;;
        [nN] ) message ERROR "Nothing more to do."; exit ;;
        * ) message ERROR "Unknown response.";;
    esac
done

message INFO "There is a preview of your paired fastq files."
COMMAND="head -n3 paired_dataset_paths.txt"
message CMD "${COMMAND}"
eval ${COMMAND}

read -p "Provide a space separated list of suffixes to remove from fastq file names to guess multiple sequencing (regex allowed, order does not matter): " regex_list
COMMAND="python3 ${PIPELINE_PREFIX}/../pyton/pair_guesser.py -e ${regex_list}"
message CMD "${COMMAND}"
eval ${COMMAND}

message INFO "Process over."