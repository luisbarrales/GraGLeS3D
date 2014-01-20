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

#Check if we should plot the network or just a grain and set PATTERN accordingly
if [ -z "$3" ]; then
	PATTERN=./Contourline_*_Timestep_
else
	if [ "$3" == "network" ]; then
		PATTERN=./Contourline_*_Timestep_
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
GNUPLOT_STRING="gnuplot -e \"set term gif animate delay 40 enhanced; set output '$1'"
while true; do
	CURRENT_PATTERN=$PATTERN$TIMESTEP
	CURRENT_PATTERN=$CURRENT_PATTERN".gnu"
	FIRST=""

	if [ "$3" == "grain" ]; then	
		if [ ! -f $CURRENT_PATTERN ]; then
			break
		fi
		if [ -z "$FIRST" ]; then 
			GNUPLOT_STRING="$GNUPLOT_STRING; plot '$CURRENT_PATTERN' w l palette title 'Timestep $TIMESTEP'"
		else
			GNUPLOT_STRING="$GNUPLOT_STRING; plot '$CURRENT_PATTERN' w l palette notitle"
		fi
	else
		if [ "$CURRENT_PATTERN" != "$(echo $CURRENT_PATTERN)" ]; then
			GNUPLOT_STRING="$GNUPLOT_STRING; plot"
			for name in $CURRENT_PATTERN; do
				if [ -z "$FIRST" ]; then
					GNUPLOT_STRING="$GNUPLOT_STRING '$name' w l palette title 'Timestep $TIMESTEP'"
					FIRST="done"
				else
					GNUPLOT_STRING="$GNUPLOT_STRING, '$name' w l palette notitle"
				fi
			done
		else
			break;
		fi
	fi
	TIMESTEP=`expr $TIMESTEP + $2`
done
echo "here"
GNUPLOT_STRING="$GNUPLOT_STRING\""
eval $GNUPLOT_STRING
echo "Done. Creted file $1."
