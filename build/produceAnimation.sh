if [ $# -lt 2 ]; then
	echo "Specify output filename and timestep size!"
	echo "Usage: produceAnimation.sh <output_filename.gif> <timestep_size> [network grain] [grain_number]"
	exit 2
fi
#Dunno how to negate this ...
if [[ $2 =~ ^-?[0-9]+$ ]]; then
	echo
else
	echo "Specify correct timestep size!"
	echo "Usage: produceAnimation.sh <output_filename.gif> <timestep_size> [network grain] [grain_number]"
	exit 2
fi
#Make a pipe because the input string will be significantly long
GNUPLOTPIPE=/tmp/gnuplotpipe
mkfifo $GNUPLOTPIPE
#Open the pipe as a file for this terminal
exec 3<> $GNUPLOTPIPE
#Check if we should plot the network or just a grain and set PATTERN accordingly
if [ -z "$3" ]; then
	PATTERN=./Network_Timestep_
else
	if [ "$3" == "network" ]; then
		PATTERN=./Network_Timestep_
	elif [ "$3" == "grain" ]; then
		if [ -z "$4" ]; then
			echo "Specify grain number!"
			echo "Usage: produceAnimation.sh <output_filename.gif> [network grain] [grain_number]"
			exit 2
		fi

		PATTERN=./Contourline_
		PATTERN=$PATTERN$4
		PATTERN=$PATTERN"_Timestep_"
	fi
fi

TIMESTEP=$2
gnuplot < $GNUPLOTPIPE &
echo "set term gif animate delay 40 enhanced" >&3
echo "set output '$1'" >&3
while true; do
	CURRENT_PATTERN=$PATTERN$TIMESTEP
	CURRENT_PATTERN=$CURRENT_PATTERN".gnu"
	FIRST=""

	if [ "$3" == "grain" ]; then	
		if [ ! -f $CURRENT_PATTERN ]; then
			break
		fi
		if [ -z "$FIRST" ]; then 
			GNUPLOT_STRING="plot '$CURRENT_PATTERN' w l palette title 'Timestep $TIMESTEP'"
		else
			GNUPLOT_STRING="plot '$CURRENT_PATTERN' w l palette notitle"
		fi
		echo "$GNUPLOT_STRING" >&3
	else
		if [ -f  $CURRENT_PATTERN ]; then			
			GNUPLOT_STRING="plot"
			GNUPLOT_STRING="$GNUPLOT_STRING '$CURRENT_PATTERN' w l palette title 'Timestep $TIMESTEP'"
			FIRST="done"
			echo "$GNUPLOT_STRING" >&3
		else
			break;
		fi
	fi
	TIMESTEP=`expr $TIMESTEP + $2`
done
GNUPLOT_STRING="$GNUPLOT_STRING\""
exec 3>&-
rm -f $GNUPLOTPIPE
echo "Done. Creted file $1."
