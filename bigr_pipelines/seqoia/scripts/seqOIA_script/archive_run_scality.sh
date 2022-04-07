#!/bin/bash
set -eo pipefail

aws_opt="--endpoint http://10.172.104.11:12290 --profile spim-preprod"
aws_prefix=""
DEBUG=${DEBUG:-1}

if [ $DEBUG -eq 1 ]; then
  aws_opt=""
  aws_prefix="echo"
fi


path=$1
if [ -z "$path" -o "$path" == "-h" -o "$path" == "--help" ]; then
  echo "Usage: $0 <path>"
  echo "Copy <path> to scality using the new tree structure"
  exit
fi

if [ ! -d "$path" ]; then
  echo "Error $path isn't a directory"
  exit 1
fi

run=$(basename $path)

json=$path/$run.json
if [ ! -e "$json" ]; then
  echo "Error can't find json $json"
  exit 1
fi

pid_crc=$(jq -r .pid_crc $json)
if [ -z "$pid_crc" ]; then
  echo "Error can't find pid_crc in $json"
  exit 1
fi
if [ $DEBUG -eq 1 ]; then
  echo $pid_crc
fi
pid_crc_2by2=$(sed -e 's/.\{2\}/&\//g' -e 's/\/$//' <<< $pid_crc)

run_name=$(sed 's/\(.*\)_\([0-9]\{2\}\)\([0-9]\{2\}\)\([0-9]\{4\}\)/\4\3\2_\1/' <<< $run)
if [ $DEBUG -eq 1 ]; then
  echo $run_name
fi

preindication=$(jq -r .preindication.catkey $json)
if [ -z "$preindication" ]; then
  echo "Error can't find preindication in $json"
  exit 1
fi
if [ $DEBUG -eq 1 ]; then
  echo $preindication
fi

# Samples
patients=$(jq -r .patients[].id_anon $json)
if [ $preindication == "p2" ]; then # cas cancer
  patients=$(jq -r ".run_id | keys | .[]" $json | xargs -i echo "${patients}_{}")
fi
if [ -z "$patients" ]; then
  echo "Error can't find any patient in $json"
  exit 1
fi

for p in ${patients}; do
  pp=$path/${p}_S*/
  extended=$(basename $pp)
  if [[ -d "$path/$extended/chr_1" ]]
  then
    for bam in $path/$extended/chr_*; do
      $aws_prefix aws $aws_opt s3 sync --exclude "*markedup*" $bam s3://spim-preprod-sample/$pid_crc_2by2/$p/$run_name/bams/;
    done
  fi
  if [ -f "$path/$extended/${extended}_identito.csv" ]; then
    $aws_prefix aws $aws_opt s3 cp $path/$extended/${extended}_identito.csv s3://spim-preprod-sample/$pid_crc_2by2/$p/$run_name/identito/
  fi
  if [ -f "$path/${extended}_identito.json" ]; then
    $aws_prefix aws $aws_opt s3 cp $path/${extended}_identito.json s3://spim-preprod-sample/$pid_crc_2by2/$p/$run_name/identito/
  fi
done

# Rest of the prescription
for f in $path/$run*.{vcf,tar.gz,csv,json}; do
  $aws_prefix aws $aws_opt s3 cp $f s3://spim-preprod-prescription/$pid_crc_2by2/$pid_crc/$run_name/
done
$aws_prefix aws $aws_opt s3 sync $path/config s3://spim-preprod-prescription/$pid_crc_2by2/$pid_crc/$run_name/
