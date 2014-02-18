FIRST=0
if [ $# -lt 1 ]; then
	echo "Timestep not specified!"
	echo "Usage: plotNetworkAtTime.sh <timestep>"
	exit 2
fi

COMMAND="gnuplot -persist -e \" plot './Network_Timestep_$1.gnu' w l palette title 'Timestep $1'\""
echo $COMMAND
eval $COMMAND
