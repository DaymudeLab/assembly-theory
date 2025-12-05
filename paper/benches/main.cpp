
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#include <csignal>              // for signal, SIGINT
#include <iostream>             // for operator<<, basic_ostream, cout, std
#include <string>               // for char_traits, allocator, operator==
#include "assemblyCalculator.h" // for assemblyCalculator
#include "help.h"               // for help
#include "ioflag.h"             // for flagParser
#include "signalHandler.h"      // for signalHandler

using namespace std;

int main(int argc, char **argv)
{
#ifdef _WIN32
    if (SetConsoleCtrlHandler((PHANDLER_ROUTINE)CtrlHandler, TRUE))
    {
    }
#else
    signal(SIGINT, signalHandler);
#endif
    if (argc > 1)
    {
        string s = argv[1];

        if (argc > 2)
            flagParser(argc, argv);

        if (s == "--help")
            help();
        else
            assemblyCalculator(s);
    }
    else
        cout << "no file selected\n";
}
