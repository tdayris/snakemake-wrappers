#!/bin/bash -e

# Authors: Thibault Dayris
# Copyright: Copyright 2021, Thibault Dayris
# Email: thibault.dayris@gustaveroussy.fr
# License: MIT

# This script is used to export functions with source

# This function only changes echo headers
# for user's sake.
function message() {
  # Define local variables
  local status=${1}         # Either INFO, CMD, ERROR or DOC
  local message="${2:-1}"   # Your message

  # Classic switch based on status
  if [ ${status} = INFO ]; then
    >&2 echo -e "\033[1;36m@INFO:\033[0m ${message}"
  elif [ ${status} = CMD ]; then
    >&2 echo -e "\033[1;32m@CMD:\033[0m ${message}"
  elif [ ${status} = ERROR ]; then
    >&2 echo -e "\033[41m@ERROR:\033[0m ${message}"
  elif [ ${status} = DOC ]; then
    >&2 echo -e "\033[0;33m@DOC:\033[0m ${message}"
  elif [ ${status} = WARNING ]; then
    >&2 echo -e "\033[1;33m@WARNING:\033[0m ${message}"
  else
    error_handling ${LINENO} 1 "Unknown message type: ${status}"
  fi
}

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
