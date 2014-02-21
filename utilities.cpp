#include "utilities.h"
#include "Settings.h"
namespace utils {
    
    
    void plotGnu(const char *fileName, const char *plotfiles, int len){
        char string_gnupl_command[len + 1000];
        strcpy( string_gnupl_command, "set term x11 persist; " );
        strcat( string_gnupl_command, "set title \"");
        strcat( string_gnupl_command, fileName);
        strcat( string_gnupl_command, "\";" );
        strcat( string_gnupl_command, "set zrange[-0.1:0.1];");
        strcat( string_gnupl_command, "set xlabel 'xrange';");
        strcat( string_gnupl_command, "set ylabel 'yrange';");
        strcat( string_gnupl_command, "set zlabel 'distance';");
        strcat( string_gnupl_command, "set contour; ");
        strcat( string_gnupl_command, "set cntrparam levels 0;");
        strcat( string_gnupl_command, "unset surface;");
		strcat( string_gnupl_command, "set view map;");
        strcat( string_gnupl_command, "splot ");
        strcat( string_gnupl_command, plotfiles);
        strcat( string_gnupl_command, "; ");
        // cout << plotfiles;
        // cout << string_gnupl_command << flush;
        FILE *pipe = popen ( "gnuplot", "w" );             // Instanz von Gnuplot
        fprintf ( pipe, "%s\n", string_gnupl_command );    // F?llen der Rohrpost zu Gnuplot
        fflush ( pipe );                                   // Plotten (flush ~ Enter)
        pclose ( pipe );
        // Beenden der Rohrpost ;) am Ende
    }

        void plotContour(const char *fileName, const char *plotfiles, int len){
        char string_gnupl_command[len + 1000];
        strcpy( string_gnupl_command, "set term x11 persist; " );
        strcat( string_gnupl_command, "set title \"");
        strcat( string_gnupl_command, fileName);
        strcat( string_gnupl_command, "\";" );
        strcat( string_gnupl_command, "set xlabel 'xrange';");
        strcat( string_gnupl_command, "set ylabel 'yrange';");
//         strcat( string_gnupl_command, "set zlabel 'distance';");
//         strcat( string_gnupl_command, "set contour; ");
//         strcat( string_gnupl_command, "set cntrparam levels 0;");
//         strcat( string_gnupl_command, "unset surface;");
// 	strcat( string_gnupl_command, "set view map;");
        strcat( string_gnupl_command, "plot ");
        strcat( string_gnupl_command, plotfiles);
        strcat( string_gnupl_command, "; ");
//         cout << plotfiles;
//         cout << string_gnupl_command << endl << endl;
	    FILE *pipe = popen ( "gnuplot", "w" );             // Instanz von Gnuplot
        fprintf ( pipe, "%s\n", string_gnupl_command );    // F?llen der Rohrpost zu Gnuplot
        fflush ( pipe );                                   // Plotten (flush ~ Enter)
        pclose ( pipe );
        // Beenden der Rohrpost ;) am Ende
    }

    
    void plotGnuPNG(const char *fileName, const char *plotfiles, int len){
        char string_gnupl_command[len + 1000];
        strcpy( string_gnupl_command, "set term x11 persist; " );
        strcat( string_gnupl_command, "set output \"");
        strcat( string_gnupl_command, fileName);
        strcat( string_gnupl_command, "\";" );
        strcat( string_gnupl_command, "set zrange[-0.1:0.1];");
        strcat( string_gnupl_command, "set xlabel 'xrange';");
        strcat( string_gnupl_command, "set ylabel 'yrange';");
        strcat( string_gnupl_command, "set zlabel 'distance';");
        strcat( string_gnupl_command, "set contour; ");
        strcat( string_gnupl_command, "set cntrparam levels 0;");
        strcat( string_gnupl_command, "unset surface;");
		strcat( string_gnupl_command, "set view map;");
        strcat( string_gnupl_command, "set terminal png;");
        strcat( string_gnupl_command, "splot ");
        strcat( string_gnupl_command, plotfiles);
        strcat( string_gnupl_command, "; ");
        // 	cout << plotfiles;
        // 	cout << string_gnupl_command << flush;
        FILE *pipe = popen ( "gnuplot", "w" );             // Instanz von Gnuplot
        fprintf ( pipe, "%s\n", string_gnupl_command );    // F?llen der Rohrpost zu Gnuplot
        fflush ( pipe );                                   // Plotten (flush ~ Enter)
        pclose ( pipe );
        // Beenden der Rohrpost ;) am Ende
    }
    
    
    void PNGtoGIF(const char *fileName) {
	
        FILE *pipe = popen ( "rm GrainNetwork.gif tempGrainNetwork.mp4; ../ffmpeg -i GrainNetwork%05d.png tempGrainNetwork.mp4; ../ffmpeg -i tempGrainNetwork.mp4 -pix_fmt rgb24 -s 640x480 GrainNetwork.gif; rm tempGrainNetwork.mp4;", "r" );
        fflush ( pipe );                                   // Plotten (flush ~ Enter)
        pclose ( pipe );
    }
    
    void CreateMakeGif(char* outputFile) {
	ofstream makeGif;
	stringstream plotFile;
	makeGif.open("makeGif.gnu");
	makeGif << "set term gif animate delay 40 enhanced;"<<endl;
	makeGif << "set output \""<< outputFile << "\"; " << endl;
	for ( int i =0; i<= Settings::NumberOfParticles; i+= Settings::AnalysysTimestep){
	    plotFile.str("");
	    plotFile<< "EnergyDistribution_T" << i << ".gnu";
	      // TODO Check if file exists? 
	    makeGif << "load \"" << plotFile.str() << "\""<< endl;
	}
	makeGif.close();
	plotFile.str("");
	plotFile << "load \"makeGif.gnu\"";
	FILE *pipe = popen ("gnuplot", "w" );
	fprintf ( pipe,"%s\n",plotFile.str().c_str()); 
	fflush ( pipe );                                   // Plotten (flush ~ Enter)
	pclose ( pipe );
    }
    
    
};
