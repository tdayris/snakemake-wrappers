#!/bin/bash

# This function will take error messages and exit the program
function error_handling() {
  # Geathering input parameter (message, third parameter is optionnal)
  echo -ne "\n"
  local parent_lineno="$1"
  local code="$2"
  local message="${3:-1}"

  # Checking the presence or absence of message
  if [[ -n "$message" ]] ; then
    # Case message is present
    message ERROR "Error on or near line ${parent_lineno}:\n ${message}"
    message ERROR "Exiting with status ${code}"
  else
    # Case message is not present
    message ERROR "Error on or near line ${parent_lineno}"
    message ERROR "Exiting with status ${code}"
  fi

  # Exiting with given error code
  exit "${code}"
}

# This function only changes echo headers
# for user's sake.
function message() {
  # Define local variables
  local status=${1}         # Either INFO, CMD, ERROR or DOC
  local message="${2:-1}"   # Your message

  # Classic switch based on status
  if [ ${status} = INFO ]; then
    echo -e "\033[1;36m@INFO:\033[0m ${message}"
  elif [ ${status} = CMD ]; then
    echo -e "\033[1;32m@CMD:\033[0m ${message}"
  elif [ ${status} = ERROR ]; then
    echo -e "\033[41m@ERROR:\033[0m ${message}"
  elif [ ${status} = DOC ]; then
    echo -e "\033[0;33m@DOC:\033[0m ${message}"
  else
    error_handling ${LINENO} 1 "Unknown message type"
  fi
}

function task(){
  sleep 0.5; echo ${1};
}

function help_message() {
  message DOC "Hi, thanks for using me as your script for parallel N-process "
  message DOC "job batch launcher. I'm very proud to be your script today,"
  message DOC "and I hope you'll enjoy working with me."
  echo ""
  message DOC "Every time you'll see a line starting with '@', it will be because I speak."
  message DOC "In fact, I always start my speech with '\033[0;33m@DOC:\033[0m' when i't about my functions,"
  message DOC "'\033[1;36m@INFO:\033[0m' when it's about my intentions, "
  message DOC "with '\033[1;32m@CMD:\033[0m' I give you all the command lines I use,"
  message DOC "and with '\033[41m@ERROR:\033[0m', I tell you when things go wrong."
  echo ""
  message DOC "I understand very fiew things, and here they are:"
  message DOC "-p | --process     The number of process per batch."
  message DOC "-f | --function    The functon i'll run in each process."
  message DOC "-d | --demo        No argument. I'll show you how I work."
  message DOC "-h | --help        Print this help message, then exit."
  echo ""
  message DOC "A typical command line would be:"
  message DOC "bash run_batch.sh -d -p 5 {a..k}"
  message DOC "This will print letters, from 'a' to 'k', batch after "
  message DOC "batch, letters five by five."
  echo ""
  message DOC "Another example:"
  message DOC 'my_samtools() {samtools index "${1}"}'
  message DOC "export -f my_samtools"
  message DOC "bash run_batch.sh -p 3 -f my_samtools *.bam"
  message DOC "This will index your bam files, batch after batch, "
  message DOC "three by three."
  exit 0
}

# Global variables
# Maximum number of process per batch
PROCESS=1
# Name of the function that is to be launched
FUNCTION=""
# Demonstration without any work
DEMO=false
# Positional arguments
POSITIONAL=()

while [ "$#" -gt 0 ]; do
  # Parse command line
  case "${1}" in
    -p|--process) PROCESS="${2}"; shift 2;;
    -f|--function) FUNCTION="${2}"; shift 2;;
    -d|--demo) DEMO=true; shift;;
    -h|--help) help_message; exit 1;;
    -*) error_handling ${LINENO} 1 "Unknown arguments ${1}"; exit 1;;
    *) POSITIONAL+=("${1}"); shift;;
  esac
done

if [ ${DEMO} = true ]; then
  FUNCTION=task
fi

if [ ${FUNCTION} == "" ]; then
  error_handling ${LINENO} 2 "Unknown function"
fi

i=0
(
for thing in "${POSITIONAL[@]}"; do
   ((i=i%PROCESS)); ((i++==0)) && wait
   message INFO "${FUNCTION} ${thing}"
   ${FUNCTION} "${thing}" &
done; wait
)
