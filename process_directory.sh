#!/bin/sh
extract_files(){
    rep="$1"
    repname="${1##*/}"
    echo "input_directory: ${repname}"
    if [ ! -d $repname ]; then
        mkdir $repname
        mkdir "$repname/waveform"
        mkdir "$repname/measurement"
        mkdir "$repname/cut"
        mkdir "$repname/calibration"
    fi 
    modifier="$2"
    for item in $rep/*; do
        echo "current_focus: ${item}"
        regex="Run_([^_]+_)+Data"
        if [ -d "$item" ]; then
            temp="${item##*/}"
            name="${temp#Run_}.root"
            echo "processing_directory: ${name}"
            for subitem in $item/*; do
                temp="${subitem##*/}"
                if [[ $temp =~ $regex ]]; then
                    temp="${BASH_REMATCH%%_Data}"
                    subname="${temp#Run_}.root"
                    echo "processing_file: ${subname}"
                    echo "surfer_girl: "
                   ./surfer_girl -in ${subitem} -out "${repname}/waveform/${name}"
                fi
            done
            echo "tree_flip: "
           ./tree_flip -in "${repname}/waveform/${name}" -out "${repname}/measurement/${name}" -mod "$modifier"
        else
             if [[ $item =~ $regex ]]; then
                temp="${BASH_REMATCH%_Data}"
                name="${temp#Run_}.root"
                echo "processing: ${name}"
                echo "surfer_girl: "
                ./surfer_girl -in ${item} -out "${repname}/waveform/${name}"
                echo "tree_flip: "
                ./tree_flip -in "${repname}/waveform/${name}" -out "${repname}/measurement/${name}" -mod "$modifier"
            fi
        fi
    done
}

if [ -z "$2" ]; then echo "second argument should contain the list of modifiers to apply: {[0-8]:[abcrt]...}..."; fi
rep="${1:-.}"
if [ "${rep: -1}" == "/" ]; then rep="${rep%/}"; fi
[ -d "$rep" ] && [ -n "$2" ] && extract_files ${rep} ${2}
