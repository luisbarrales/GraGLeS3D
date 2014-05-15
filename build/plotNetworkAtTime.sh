FIRST=0
if [ $# -lt 1 ]; then
	echo "Timestep not specified!"
	echo "Usage: plotNetworkAtTime.sh <timestep>"
	exit 2
fi

COMMAND="gnuplot -persist -e \"set palette rgbformulae 33,13,10; set cbrange[0:0.6]; set cbtics 0.1; set yrange [:] reverse; plot './Network_Timestep_$1.gnu' w l palette title 'Timestep $1'\""
echo $COMMAND
eval $COMMAND
