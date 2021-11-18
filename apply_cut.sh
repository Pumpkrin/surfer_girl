#!/bin/sh
apply_cut(){
    file="$1"
    temp="${1##*/}"
    filename="${temp%.cut}.root"
    echo $filename
    
    root_directory="${1%%/cut/*}"
    modifier="$2"
    echo "processing: ${name}"
    echo "tree_grubbing: "
    ./tree_grubbing -cut $file 
    echo "tree_flip: "
    ./tree_flip -in "${root_directory}/waveform/${filename}" -out "${root_directory}/measurement/${filename}" -mod "$modifier"
}

if [ -z "$2" ]; then echo "second argument should contain the list of modifiers to apply: {[0-8]:[abcrt]...}..."; fi
[ -f "$1" ] && [ -n "$2" ] && apply_cut ${1} ${2}
