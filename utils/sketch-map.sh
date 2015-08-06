#!/bin/bash

# modify this variable so that it points to the path where the sketch-map
# executables are stored
SMAPROOT=$HOME/bin/

# makes sure we have a valid dimred path
if [ ! -e "$SMAPROOT/dimred" ]; then
   read -p "Couldn't find dimred. Please enter the path where dimred resides " SMAPROOT
fi
SMAP="$SMAPROOT/dimred"

read -p "Please enter the dimensionality of input data " HD
LD=2      #hardcoded, dimred only works with d=2 right now.

read -p "Are we reading the similarity matrix? " SIM
if [ $SIM = "y" ]; then SIM=" -similarity "; else SIM=""; fi

read -p "Are points weighted [y/n]? " DOW
if [ $DOW = "y" ]; then DOW=" -w "; else DOW=""; fi

read -p "Should I use dot-product distance [y/n]? " DOT
if [ $DOT = "y" ]; then DOT=" -dot "; else DOT=""; fi


PI=0
read -p "Please enter the periodicity of input data [0 if non-periodic] " PI
if [ -z $PI -o $PI = "0" ]; then PI=""; else  PI=" -pi $PI"; fi

read -p "Please enter the input data file name " FILEHD
read -p "Please enter the output data prefix " FILELD

read -p "Please enter high dimension sigma, a, b [e.g. 6.0 2 6 ] " SIGMAHD AHD BHD
read -p "Please enter low  dimension sigma, a, b [e.g. 6.0 2 6 ] " SIGMALD ALD BLD

rm log
echo "Now running a preliminary iterative metric MDS and sketch-map."
if [ ! -e $FILELD.imds ]; then
   echo "$SMAP -vv -D $HD -d $LD $PI $DOW $DOT $SIM -center -preopt 100"
   grep -v \#  $FILEHD | $SMAP -vv -D $HD -d $LD $PI $DOW $DOT $SIM -center -preopt 100 > $FILELD.imds 2>>log
fi
grep -v "#" $FILELD.imds | awk '{print $1, $2}' > tmp
if [ ! -e $FILELD.ismap ]; then
   echo "$SMAP -vv -D $HD -d $LD $PI $DOW $DOT $SIM -center -preopt 100 -fun-hd $SIGMAHD,$AHD,$BHD -fun-ld $SIGMALD,$ALD,$BLD -init tmp"
   grep -v \#  $FILEHD | $SMAP -vv -D $HD -d $LD $PI $DOW $DOT $SIM -center -preopt 100 -fun-hd $SIGMAHD,$AHD,$BHD -fun-ld $SIGMALD,$ALD,$BLD -init tmp > $FILELD.ismap 2>> log
fi

GW=`awk 'BEGIN{maxr=0} !/#/{r=sqrt($1^2+$2^2); if (r>maxr) maxr=r} END{print maxr*1.2}' $FILELD.imds`;
NERR=`awk '/Error/{print $(NF)}'  $FILELD.imds | tail -n 1`
SMERR=`awk '/Error/{print $(NF)}'  $FILELD.ismap | tail -n 1`

IMIX=1.0
MAXITER=10
for ((ITER=1; ITER<=$MAXITER; ITER++)); do
   MDERR=$NERR
   echo "Mixing in $IMIX"
   if [ ! -e $FILELD.gmds_$ITER ]; then
      echo "Now running $SMAP -vv -D $HD -d $LD $PI $DOW $DOT $SIM -center -preopt 50 -grid $GW,21,201 -fun-hd $SIGMAHD,$AHD,$BHD -fun-ld $SIGMALD,$ALD,$BLD -init tmp -gopt 3 -imix $IMIX < $FILEHD > $FILELD.gmds_$ITER 2>>log"
      grep -v \#  $FILEHD | $SMAP -vv -D $HD -d $LD $PI $DOW $DOT $SIM -center -preopt 50 -grid $GW,21,201 -fun-hd $SIGMAHD,$AHD,$BHD -fun-ld $SIGMALD,$ALD,$BLD -init tmp -gopt 3 -imix $IMIX > $FILELD.gmds_$ITER 2>>log
   fi
   grep -v "#" $FILELD.gmds_$ITER | awk '{print $1, $2}' > tmp
   GW=`awk 'BEGIN{maxr=0} !/#/{r=sqrt($1*$1+$2*$2); if (r>maxr) maxr=r} END{print maxr*1.2}' $FILELD.gmds_$ITER`;
   NERR=`awk '/Error/{print $(NF)}' $FILELD.gmds_$ITER  | tail -n 1 `
   echo "Residual error is $NERR"
   IMIX=`echo "$IMIX  $SMERR  $NERR" | awk '{new=$2/($2+$3); if (new<0.1) new=0.1; if (new>0.5) new=0.5; print new*$1 }'`
   echo "DEBUG  $MDERR $NERR"
   if [ ` echo $MDERR $NERR | awk -v i=$ITER '{ if (i>1 && (($1-$2)/$2)*(($1-$2)/$2)<1e-4) print "done"; else print "nope";}' ` = "done" ]; then ((ITER++)); break; fi;
done

echo "Doing final fit"
((ITER--))
grep -v "#" $FILELD.gmds_$ITER | awk '{print $1, $2}' > tmp
grep -v \#  $FILEHD | $SMAP -vv -D $HD -d $LD $PI $DOW $DOT $SIM -center -preopt 100 -grid $GW,21,201 -fun-hd $SIGMAHD,$AHD,$BHD -fun-ld $SIGMALD,$ALD,$BLD -init tmp -gopt 10 > $FILELD.gmds 2>>log

