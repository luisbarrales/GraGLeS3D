if [ $# -lt 1 ]; then
	echo "Specify output filename and timestep size!"
	echo "Usage: producePNGs.sh <timestep_size>"
	exit 2
fi
#Dunno how to negate this ...
if [[ $1 =~ ^-?[0-9]+$ ]]; then
	echo
else
	echo "Specify correct timestep size!"
	echo "Usage: produceAnimation.sh <output_filename.gif> <timestep_size>"
	exit 2
fi
PATTERN=./Network_Timestep_
TIMESTEP=0
while true; do
	CURRENT_FILENAME=$PATTERN`printf "%06d" $TIMESTEP`.png
	echo -e -n "Processing timestemp $TIMESTEP\r"
	CURRENT_PATTERN=$PATTERN$TIMESTEP
	CURRENT_PATTERN=$CURRENT_PATTERN".gnu"
	FIRST=""
	if [ -f  $CURRENT_PATTERN ]; then			
		GNUPLOT_STRING="gnuplot -e \"set term png giant size 1024,768; set output '$CURRENT_FILENAME'; plot"
		GNUPLOT_STRING="$GNUPLOT_STRING '$CURRENT_PATTERN' w l palette title 'Timestep $TIMESTEP' \""
		FIRST="done"
	eval $GNUPLOT_STRING
	else
		break;
	fi
	TIMESTEP=`expr $TIMESTEP + $1`
done
echo ""
echo "Done creating intermidiate files."
echo "Creating gif with resolution 1024 , 768"
convert -delay 40 -loop 0 Network_Timestep_*.png NetworkEvolution.gif
echo "Done creating gif. Cleaning up.."
rm *.png
echo "Done."

