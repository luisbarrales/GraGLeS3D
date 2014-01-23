PATTERN=./Contourline_*_Timestep_$1.gnu
COMMAND=plot
GNUPLOT="gnuplot -persist -e "
FIRST=0
if [ $# -lt 1 ]; then
	echo "Timestep not specified!"
	echo "Usage: plotNetworkAtTime.sh <timestep>"
	exit 2
fi

for name in $PATTERN; do
	if [ "$FIRST" -eq  0 ]; then
		if [ "$name" = "$PATTERN" ]; then
			echo "Files for the selected timestep do not exist !"
			exit 2
		fi
		COMMAND="$COMMAND '$name' w l palette title 'Timestep $TIMESTEP'"
		FIRST=1
	else
		COMMAND="$COMMAND, '$name' w l palette notitle"
	fi
done
GNUPLOT="$GNUPLOT \"$COMMAND\""
eval $GNUPLOT
