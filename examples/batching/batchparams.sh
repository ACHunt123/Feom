#!/bin/bash

### patameter files for the batch run
script="~/software/phd/Feom/feom.py"


# Path to your params file
PARAMS_FILE="batchparams.json"

#copy the params file to plotparams.json5 for plotting later
# the json5 package allows comments in the json file (for easier plotting)
OUTFILE="plotparams.json5"

awk '
/:/ {
  # Find the first colon in the line
  colon_pos = index($0, ":")
  before = substr($0, 1, colon_pos)
  after = substr($0, colon_pos + 1)

  # Clean up comment: trim leading spaces
  comment = after
  sub(/^[ \t]+/, "", comment)

  # Add comment at the end
  print before after " // " comment
  next
}
{ print }
' "$PARAMS_FILE" > "$OUTFILE"


# Read arrays using jq
Ks=($(jq '.Ks[]' "$PARAMS_FILE"))
Ls=($(jq '.Ls[]' "$PARAMS_FILE"))
dts=($(jq '.dts[]' "$PARAMS_FILE"))
etas=($(jq '.etas[]' "$PARAMS_FILE"))
betas=($(jq '.betas[]' "$PARAMS_FILE"))


# export MPLCONFIGDIR=/scratch2/ach221/.config/matplotlib
# HOME=/scratch2/ach221
# script="$HOME/software/phd/Feom/feom.py"


# hardcoded parameters using jq
HRDparams=$(jq -r '.HRDparams | to_entries | map("--\(.key) \(.value)") | join(" ")' "$PARAMS_FILE")


#number of parrallel processes allowed 
nparallel=30

#whether to quiet the output or not
quiet=''
quiet='> /dev/null 2>&1'