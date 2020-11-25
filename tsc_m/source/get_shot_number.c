#include <stdio.h>
#include <string.h>

int get_shot_number_(char *transp_number)
{
int shot_number;

    shot_number = atoi(transp_number);
/*
    printf( "transp_number=%s  shot_number=%d\n", 
             transp_number,    shot_number ); 
*/
    return shot_number;
}
