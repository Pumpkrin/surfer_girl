#!/bin/sh
extract_file(){
    file="$1"
    filename="${1##*/}"
    directory="temp"
    if [ ! -d "${directory}" ]; then
        mkdir "${directory}"
        mkdir "${directory}/waveform"
        mkdir "${directory}/measurement"
        mkdir "${directory}/cut"
        mkdir "${directory}/calibration"
    fi 
    modifier="$2"
    regex="Run_([^_]+_)+Data"
    if [[ $file =~ $regex ]]; then
        temp="${BASH_REMATCH%_Data}"
        name="${temp#Run_}.root"
        echo "processing: ${name}"
        echo "surfer_girl: "
        ./surfer_girl -in ${file} -out "${directory}/waveform/${name}"
        echo "tree_flip: "
        ./tree_flip -in "${directory}/waveform/${name}" -out "${directory}/measurement/${name}" -mod "$modifier"
    fi
}

if [ -z "$2" ]; then echo "second argument should contain the list of modifiers to apply: {[0-8]:[abcrt]...}..."; fi
[ -f "$1" ] && [ -n "$2" ] && extract_file ${1} ${2}
