#/bin/bash

AGGREGATE="/home/jpablo/mis_bin"

#DEFAULTS
inp="aggregate.inp"
top="topol.top"
xtc="md.xtc"
out="output.dat"
#
schm="original"
scrn="complete"
conf="body"
doscrn="NO"
doconf="NO"
dopim="NO"
add="[none]"
#
size=10
#
bin="aggregate"
np="1"
help="NO"

#TAKE INITIAL COMMAND BEFORE PROCESSING IT
initial_command=$@

# READING INPUT DATA FROM COMMAND LINE
while test "x${1}" != x ; do
    case ${1} in
     -f         ) shift; inp="${1}"   ;;
     -p         ) shift; top="${1}"   ;;
     -t         ) shift; xtc="${1}"   ;;
     -o         ) shift; out="${1}"   ;;
     -m         ) shift; size="${1}"  ;;
     -b         ) shift; bin="${1}"   ;;
     -np        ) shift; np="${1}"    ;;
     -r         ) shift; schm="${1}"  ;;
     -s         ) shift; scrn="${1}" ; doscrn="YES"  ;;
     -c         ) shift; conf="${1}" ; doconf="YES"  ;;
     -a         ) shift; add="${1}"   ;;
     -life      ) dolife="YES"        ;;
     -nolife    ) dolife="NO"         ;;
     -scrn      ) doscrn="YES"        ;;
     -noscrn    ) doscrn="NO"         ;;
     -conf      ) doconf="YES"        ;;
     -noconf    ) doconf="NO"         ;;
     -pim       ) dopim="YES"         ;;
     -nopim     ) dopim="NO"          ;;
     -h         ) help="YES"          ;;
     *          ) bad="$bad ${1}" ; error="YES"
    esac
    shift
done

if [ "$bad" != "" ]; then
    echo
    echo ERROR: Unknown statements from command line: $bad
    echo
    exit 1
fi

if [ "${inp: -4}" != ".inp" ] ; then
  inp=${inp}.inp
  echo "NOTE: Added default extension to input file name: ${inp}"
  echo
fi

if [ "${top: -4}" != ".top" ] ; then
  top=${top}.top
  echo "NOTE: Added default extension to topology file name: ${top}"
  echo
fi

if [ "${xtc: -4}" != ".xtc" ] ; then
  xtc=${xtc}.xtc
  echo "NOTE: Added default extension to trajectory file name: ${xtc}"
  echo
fi
if [ "${out: -4}" != ".dat" ] ; then
  out=${out}.dat
  echo "NOTE: Added default extension to output file name: ${out}"
  echo
fi

name="${out:0:-4}"

#PRINT PARAMETERS IN USE
cat <<EOF

=========================================================
INPUT DATA
   -f            input file name                     ${inp}
   -c            topology file name                  ${top}
   -t            trajectory file name                ${xtc}
   -o            output file name                    ${out}

ALGORITHM OPTIONS
   -m            maximum aggregate size             ${size}
   -r            geometric criteria scheme          ${schm}
   -[no]scrn     activate screening algorithm       ${doscrn}
   -[no]life     activate lifetime algorithm        ${dolife}
   -[no]conf     activate conformational analysis   ${doconf}
   -[no]pim      pairwise interaction matrix        ${dopim}

   -s            screening algorithm                ${scrn}
   -c            conformational analysis algorithm  ${conf}

ADDITIONAL OPTIONS
   -b            binary file name                   ${bin}
   -np           number of procesors per            ${np}
                  node to use
   -a            additional keyword(s)              ${add}

OPTIONS
   -h            print help and quit                ${help}
=========================================================

EOF

#EXITING AFTER HELP
if [ "$help" == "YES" ]; then
  exit 0
fi

# GENERAL SETTINGS AND FATAL ERRORS CHECK

if [ "${dopim}" == "YES" ]; then
  pim="-pim"
elif [ "${dopi,}" == "NO" ]; then
  pim="-nopim"
fi

if [ "${dolife}" == "YES" ]; then
  life="-life"
elif [ "${dolife}" == "NO" ]; then
  life="-nolife"
fi

if [ "${doscrn}" == "YES" ]; then
  scrn="-s ${scrn}"
elif [ "${doscrn}" == "NO" ]; then
  scrn="-noscrn"
fi

if [ "${doconf}" == "YES" ]; then
  conf="-c ${conf}"
elif [ "${doconf}" == "NO" ]; then
  conf="-noconf"
fi

if [ "${add}" == "[none]" ]; then
  add=""
fi 

aggregate=${bin}
#aggregate=${AGGREGATE}/${bin}
command -v $aggregate >/dev/null || { echo "ERROR: Cannot find aggregate binary: ${aggregate}"; exit 1; }

test -f $inp || echo "ERROR: Input file \"$inp\" not found in current folder"
test -f $inp || exit 1

test -f $top || echo "ERROR: Topology file \"$top\" not found in current folder"
test -f $top || exit 1

test -f $xtc || echo "ERROR: Trajectory file \"$xtc\" not found in current folder"
test -f $xtc || exit 1

# Preparing the script to send the calculation
cat > aggregate_${name} <<!
#!/bin/bash

export OMP_NUM_THREADS=${np}

${aggregate} -np ${np} -f ${inp} -p ${top} -t ${xtc} -o ${out} -m ${size} -r ${schm} ${pim} ${life} ${scrn} ${conf} > ${out//.dat/.out}

!

chmod u+x qsub_${name}_aggregate

qsub qsub_${name}_aggregate -np ${np}
