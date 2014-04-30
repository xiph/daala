#!/bin/sh
#Usage: make_testcase.sh <output directory> <directory of reference images> \
# <directory of images coded with condition 1> \
# <directory of images coded with condition 2> ...
#The directory for each condition must contain .png images with the same names
# as those in the directory of reference images.
#All paths must be on the same volume (to allow hard links to be created in the
# output directory).
#Creates an index.html in the output directory that shows the images with the
# order of the conditions randomized for each image, and an <output>-key.txt
# file with the key to unscramble the results.
output="${1%%/}"
shift
references="${1%%/}"
shift
nconditions=$#
declare -a conditions
for ((i=0;i<nconditions;i++)) ; do
  conditions[i]="${1%%/}"
  shift
done
echo "nconditions: $nconditions"

function randn () {
  local n=$1
  local m=$((32768-(32768%n)))
  local ret=$m
  while [[ "$ret" -ge "$m" ]] ; do
    ret=${RANDOM}
  done
  echo $((ret%n))
}

function padidx () {
  local i="${1}"
  local idx="${i}"
  [[ "${i}" -lt 10000 ]] && idx="0${idx}"
  [[ "${i}" -lt 1000 ]] && idx="0${idx}"
  [[ "${i}" -lt 100 ]] && idx="0${idx}"
  [[ "${i}" -lt 10 ]] && idx="0${idx}"
  echo "${idx}"
}

declare -a ref_imgs
declare -a test_imgs
declare -a test_cond
nrefs=0
for f in "${references}"/*.png ; do
  base="${f##*/}"
  ref_imgs[nrefs]="${base}"
  for ((i=0;i<nconditions;i++)) ; do
    if [[ ! -f "${conditions[i]}/${base}" ]] ; then
      echo "Missing ${conditions[i]}/${base}"
      break
    fi
  done
  [[ "$i" -ge "${nconditions}" ]] && nrefs=$((nrefs+1))
done
echo "nrefs: $nrefs"

mkdir -p "${output}"

index="${output}/index.html"
key="${output}-key.txt"
rm "${key}" > /dev/null 2>&1

cp -a "header.html" "${index}"
echo '<DIV CLASS="hidden"><SCRIPT LANGUAGE="javascript" TYPE="text/javascript">' >> "${index}"
echo '<!--//--><![CDATA[//><!--' >> "${index}"
echo 'function preloader() {' >> "${index}"
for ((i=0;i<nrefs;i++)) ; do
  idx=`padidx "${i}"`
  echo '  preloadImage('\'"${idx}ref.png"\'');' >> "${index}"
  for ((j=0;j<nconditions;j++)) ; do
    echo '  preloadImage('\'"${idx}-${j}"\'');' >> "${index}"
  done
done
echo '}' >> "${index}"
echo 'window.onload=preloader;' >> "${index}"
echo '//--><!]]>' >> "${index}"
echo '</SCRIPT></DIV>' >> "${index}"
echo '<TR><TD COLSPAN="'"$((nconditions+2))"'"></TD><TD><DIV ID="img_label">&nbsp;</DIV></TD></TR>' >> "${index}"
echo '<TR><TD COLSPAN="'"$((nconditions+2))"'"></TD><TD VALIGN="top" ROWSPAN="'$((nrefs+1))'"><IMG ID="img"></TD></TR>' >> "${index}"
for ((i=0;i<nrefs;i++)) ; do
  idx="`padidx ${i}`"
  j="`randn $((nrefs-i))`"
  j="$((i+j))"
  base="${ref_imgs[j]}"
  ref_imgs[j]="${ref_imgs[i]}"
  ref_imgs[i]="${base}"
  test_imgs[0]="${codeca}/${base}"
  test_imgs[1]="${codecb}/${base}"
  cp -alf "${references}/${base}" "${output}/${idx}ref.png"
  echo -n "${idx} " >> "${key}"
  echo -n '<TR><TD>'"$i"'</TD><TD><A HREF="'"${idx}ref.png"'" TARGET="_blank" onMouseOver="setImage('\'"${idx}ref.png"\'')">Reference</A></TD>' >> "${index}"
  for ((j=0;j<nconditions;j++)) ; do test_cond[j]="${j}" ; done
  for ((j=0;j<nconditions;j++)) ; do
    k="`randn $((nconditions-j))`"
    k="$((j+k))"
    tmp="${test_cond[k]}"
    test_cond[k]="${test_cond[j]}"
    test_cond[j]="${tmp}"
    cp -alf "${conditions[test_cond[j]]}/${base}" "${output}/${idx}-${j}.png"
    echo -n "${conditions[test_cond[j]]} " >> "${key}"
    echo -n '<TD><A HREF="'"${idx}-${j}.png"'" TARGET="_blank" onMouseOver="setImage('\'"${idx}-${j}.png"\'')">Condition&nbsp;'"${j}"'</A></TD>' >> "${index}"
  done
  echo "${ref_imgs[i]}" >> "${key}"
  echo '</TR>' >> "${index}"
done
cat footer.html >> "${index}"
