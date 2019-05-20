#!/usr/bin/env bash
#
# Galaxy wrapper for mintia assembly
# 

set -e

export PATH=$PATH:$(dirname $0)

# Variables
zipin=$1
vectorSeq=$2
out_zip=$3
out_html=$4

## Separate significance and posterior arguments from arguments to mintia
until [ $# -eq 0 ]
do
  case $1 in
    zipin=*)
      zipin=${1#zipin=}
      ;;
    vectorSeq=*)
      vectorSeq=${1#vectorSeq=}
      ;;
    out_zip=*)
      out_zip=${1#out_zip=}
      ;;
    out_html=*)
      out_html=${1#out_html=}
      ;;
    *)
      if [ -z "$new_args" ]; then
        new_args=$1
      else
        new_args="$new_args $1"
      fi
      ;;
  esac

  shift
done



#TFILE="/tmp/MINTIA1.$$.tmp"
#TODO : define this temporary path in version 18
TFILE="/galaxydata/galaxy-preprod/my_workspace/MINTIA1.$$.tmp"

## Run mintia assembly
mkdir $TFILE; cd $TFILE/; chmod -R 777 $TFILE/; mkdir DATA/; chmod -R 777 DATA/; cp $zipin $TFILE/DATA.zip; unzip $TFILE/DATA.zip -d $TFILE/DATA/; 
mkdir $TFILE/RESULTS/; chmod -R 777 $TFILE/RESULTS/;
cd $TFILE/DATA/*/; ln -s /galaxydata/galaxy-preprod/my_tools/MINTIA1/mintia_assembly.pl .;
echo "perl mintia_assembly.pl --input '$TFILE'/DATA/*/* --vectorSeq '$vectorSeq' --dirOutputs '$TFILE'/RESULTS";
#TODO with version 18 : path to scripts to define with a variable
perl /galaxydata/galaxy-preprod/my_tools/MINTIA1/mintia_assembly.pl --input $TFILE/DATA/*/* --vectorSeq $vectorSeq --dirOutputs $TFILE/RESULTS  1>/dev/null;

if [ $? -ne 0 ]; then
  echo "failed: perl '$__tool_directory__/mintia_assembly.pl' --input *.* --vectorSeq $vectorSeq --dirOutputs $TFILE/RESULTS/ "
  exit 1
fi

#Prepare results files as Galaxy outputs
cd $TFILE; mkdir FINAL/; chmod -R 777 FINAL/; 

for fasta in `\ls -1 $TFILE/RESULTS/*/contigs.fasta`; do
                        rn1="$fasta"
                        rn2="" 
                        rn3="${rn1/$TFILE\/RESULTS\//$rn2}"
                        rn3="${rn3/\/contigs.fasta/$rn2}"
                        echo "rn3 : $rn3"
                        echo "cp $fasta $TFILE/FINAL/$rn3.contigs.fasta"
                        cp "$fasta" "$TFILE"/FINAL/"$rn3".contigs.fasta
done
  
echo "zip -r $TFILE/final.zip $TFILE/FINAL/'"
zip -r $TFILE/final.zip $TFILE/FINAL/ 1>/dev/null;
if [ $? -ne 0 ]; then
  echo "failed: zip -r '$TFILE'/final.zip '$TFILE'/FINAL/ "
  exit 1
fi

## Move outputs files
echo "cp -a '$TFILE'/final.zip '$out_zip'"
cp -a $TFILE/final.zip $out_zip;
echo "mv  '$TFILE'/RESULTS/mintia.html '$out_html'"
#mv to move html report in final directory in order to cleanup RESULTS , DATA and DATA.zip files
mv $TFILE/RESULTS/mintia.html $out_html;

## Cleanup
rm -rf $TFILE/DATA/;
rm -rf $TFILE/RESULTS/;
rm -rf $TFILE/DATA.zip;

