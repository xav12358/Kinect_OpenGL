#include <include/glWindow2.h>
#include <include/system.h>

#include "iostream"

using namespace std;

extern GLWindow2 gGLWin;


int main( int argc, char **argv )
{


    System *S = new System(gGLWin);
    S->Run();


    return 0;
}
