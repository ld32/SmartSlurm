#!/bin/sh

setSoftware() {
  [ -z "$1" ] && { echo "Usage: $0 <bowtie-r>". Return: bowtie-hg19; return; }
  tmp=${1##*-} 
  echo got ${1%-*}-${!tmp}
  tmp=${1%-*}-${!tmp} 
  eval "$2=$tmp"
}
export -f setSoftware

setInputs() {
  [ -z "$1" ] && { echo "Usage: $0 <read1.read2>". Return: bowtie-hg19; return; }
  tmp=""
  for i in ${1/./ }; do 
    tmp="${!i},$tmp"
  done
  eval "$2=\"$tmp\""
}
export -f setInputs

rsyncToTmp() {
  [ -z "$1" ] && { echo "Usage: $0 </tmp/indexPath> </tmp/gtfPath> </tmp/bowtieIndexPath> ... "; return; } 
  for v in $@; do
      echo Working to copy: $v, waiting lock...
      [[ $v == /tmp/smartSlurm/* ]] || continue
      ls ${v#/tmp/smartSlurm/}* >/dev/null 2>&1 || { echo Reference file or folder not exist: ${v#/tmp/smartSlurm}; continue; } 
      lockFile=/tmp/${v//\//-}
      
      while ! ( set -o noclobber; echo "$$" > "$lockFile") 2> /dev/null; do
        echo waiting for lock file: $lockFile
        sleep 30
      done
      echo Got lock: $lockFile. Copying data to: $v
      trap 'rm -f "$lockFile"; exit $?' INT TERM EXIT
      #echo checking:
      mkdir -p ${v%/*}
      rsync  -aL ${v#/tmp/smartSlurm/}* ${v%/*}/
      chmod -R 777 /tmp/smartSlurm/
      find  ${v%}* -type f -exec chmod -x '{}' \;
      
      echo Copying is done for $v
      rm -f "$lockFile"
      trap - INT TERM EXIT
  done
}
export -f rsyncToTmp

setPath() {
 
  [ -z "$1" ] && { echo "Usage: $0 index.gtf.bowtieIndex"; return; } 
  for v in ${1//./ }; do
      #echo Reset path for: $v
      vx=${!v}
      [[ $vx == /tmp/* ]] && continue
      #vx=`realpath $vx` 
      eval $v=/tmp/smartSlurm/$vx
  done
  #echo new path: $gtf, $index
}
export -f setPath

setPathBack() {
 
  [ -z "$1" ] && { echo "Usage: $0 index.gtf.bowtieIndex"; return; } 
  for v in ${1//./ }; do
      #echo ResetBack path for: $v
      vx=${!v}
      [[ $vx != /tmp/* ]] && continue
      #vx=`realpath $vx` 
      eval $v=${vx#/tmp/smartSlurm/}
  done
  #echo new path: $gtf, $index
}
export -f setPathBack
