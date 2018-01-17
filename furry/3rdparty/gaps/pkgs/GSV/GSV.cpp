/* Source file for GSV module */



/* Include files */

#include "GSV.h"



/* Private variables */

static int GSV_active_count = 0;



int R3InitShapes(void)
{
    // Check whether are already initialized 
    if ((GSV_active_count++) > 0) return TRUE;

    // Initialize dependencies
    if (!R3InitShapes()) return FALSE;

    // return OK status 
    return TRUE;
}



void R3StopShapes(void)
{
    // Check whether have been initialized 
    if ((--GSV_active_count) > 0) return;

    // Stop dependencies
    R3StopShapes();
}




