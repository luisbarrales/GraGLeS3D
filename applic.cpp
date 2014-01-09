#include "applic.h"
#include <stdio.h>

void exitus (const char *s)
{
        printf("%s\n",s);
		//MPI_Abort(MPI_COMM_WORLD,1);
        //exit(1);
}
