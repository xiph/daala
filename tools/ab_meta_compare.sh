#!/bin/bash
set -e

# script to generate one image or video file based on the given input
# bpp, filename, and time. It is stored in `runs/$time/$bpp`

bpp="$1"
time="$2"
workset="$3" # subset1 or video-1-short for example
file_path="$4" # is a `*.y4m` file

base_path="/home/ec2-user"
daala_path="$base_path/daala"
libvpx_path="$base_path/libvpx"
x264_path="$base_path/x264"
x265_path="$base_path/x265"
tools_path="$daala_path/tools"
output_path="$base_path/runs/$time/$workset/bpp_$bpp"

options="-b $bpp
         -c daala
         -l $libvpx_path
         -x $x264_path
         -X $x265_path
         -d $daala_path
         -a $tools_path"

echo "Encoding $file_path with a bpp of $bpp"

mkdir --parents "$output_path"
cd "$output_path"

"$tools_path/ab_compare.sh" $options "$file_path"

basename=$(basename "$file_path")

# Remove any generated .y4m files from $output_path
rm *.y4m

# Bash parameter expansion: `${foo%.*}` searches `foo` from the back looking to match `.*`
# (find something starting with a period). If it finds it, it deletes it. A `%%` would delete
# from the first period to the end.
# This returns the filename without the `.y4m` extension.
file=${basename%.*}

if [[ $workset != *"video"* ]]
then
    # Remove newly generated files ending with .ogv
    rm "$file"*.ogv

    ext="png"
else
    ext="ogv"
fi

# ab_compare renames `a.y4m` to something like `a-70.ogv.png` so
# get `a` off the original filename and append `.png` to it.
# This might have issues if the same file was output with different
# bpp values in this folder. However, this script puts different bpp
# values in different folders so it shouldn't have any issues.
mv "$file"*".$ext" "$file.$ext"


cd "$base_path"
