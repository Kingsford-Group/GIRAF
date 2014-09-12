#include <string>
#include <vector>
#include <iostream>
#include "timer.h"

#define PROG_NAME "extract_reassortments"

using namespace std;

vector<string> main_extract_reassortments(int, char**);

int
main(int argc, char *argv[])
{
    Timer T;
    cout << PROG_NAME << ": " << T.start();

    main_extract_reassortments(argc, argv);

    cout << PROG_NAME << ": " << T.stop();
    cout << PROG_NAME << ": " << T.report() << endl;
    return 0;
}
