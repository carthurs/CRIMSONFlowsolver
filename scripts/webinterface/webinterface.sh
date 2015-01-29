#!/bin/bash
SCRIPTLOCATIONDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

if test -z "$1"; then
 echo "Usage: webinterface.sh <number of surfaces in PressHist.dat>"
 exit 1
fi

surfacesToPlotUpperBound=$(($1+1))

write_rcr_data.gpi PressHist.dat
write_rcr_data.gpi FlowHist.dat
gnuplot -e "PressFlowHistSurfaces='$surfacesToPlotUpperBound'" $SCRIPTLOCATIONDIR/makehtmlfigs.gp

rm results.html
echo "<html>" >> results.html
echo " <frameset cols=\"20%,80%\">" >> results.html
echo "  <frame src=\"bar.html\">" >> results.html
echo "  <frame src=\"\" name=\"viewer\">" >> results.html
echo " </frameset>" >> results.html
echo "</html>" >> results.html

rm bar.html
echo "<html>" >> bar.html
echo " <a href=\"QMV.svg\" target=\"viewer\">QMV</a>" >> bar.html
echo "<br>" >> bar.html
echo " <a href=\"PLV.svg\" target=\"viewer\">PLV</a>" >> bar.html
echo "<br>" >> bar.html
echo " <a href=\"VLV.svg\" target=\"viewer\">VLV</a>" >> bar.html
echo "<br>" >> bar.html
echo " <a href=\"ELV.svg\" target=\"viewer\">ELV</a>" >> bar.html
echo "<br>" >> bar.html
echo " <a href=\"Paorta.svg\" target=\"viewer\">Paorta</a>" >> bar.html
echo "<br>" >> bar.html
echo " <a href=\"Qaorta.svg\" target=\"viewer\">Qaorta</a>" >> bar.html
echo "<br>" >> bar.html
echo " <a href=\"Pstab.svg\" target=\"viewer\">Pstab</a>" >> bar.html
echo "<br>" >> bar.html
echo " <a href=\"sysAhist.svg\" target=\"viewer\">sysAhist</a>" >> bar.html
echo "<br>" >> bar.html
echo " <a href=\"MVO2_history.svg\" target=\"viewer\">MVO2_history</a>" >> bar.html
echo "<br>" >> bar.html
echo " <a href=\"r_d_history1.svg\" target=\"viewer\">r_d_history1</a>" >> bar.html
echo "<br>" >> bar.html
echo " <a href=\"r_p_history1.svg\" target=\"viewer\">r_p_history1</a>" >> bar.html
echo "<br>" >> bar.html
echo " <a href=\"meanO2Discrepancy_history1.svg\" target=\"viewer\">meanO2Discrepancy_history1</a>" >> bar.html
echo "<br>" >> bar.html
for i in `seq 1 $1`; do
 echo " <a href=\"outletPressure$i.svg\" target=\"viewer\">outletPressure$i</a>" >> bar.html
 echo "<br>" >> bar.html
done;
for i in `seq 1 $1`; do
 echo " <a href=\"outletFlow$i.svg\" target=\"viewer\">outletFlow$i</a>" >> bar.html
 echo "<br>" >> bar.html
done;

firefox results.html
