#include <iostream>
#include "optparser.h"
#include "util.h"

using namespace std;

OptParser::OptParser (int argc, char** argv, const char* progname, const char* brief)
  : argvec(argc),
    prog(progname),
    unlimitImplicitSwitches(false)
{
  for (int n = 0; n < argc; ++n)
    argvec[n] = argv[n];
  argvec.pop_front();  // program path
  briefText = "Usage: " + prog + ' ' + brief + '\n';
  text = briefText + "\n";
}

string OptParser::getCommand (const char* error) {
  if (argvec.empty()) {
    if (error)
      cerr << error << endl;
    else
      cerr << briefText;
    exit (EXIT_FAILURE);
  }
  const string command (argvec[0]);
  argvec.pop_front();
  return command;
}

bool OptParser::parseUnknown() {
  if (argvec.size()) {
    string arg (argvec[0]);
    if (arg == "-abort") {
      // test stack trace
      Abort ("abort triggered");

    } else {
      if (arg[0] == '-' || implicitSwitches.empty()) {
	cerr << text << "Unknown option: " << arg << endl;
	cerr << "Error parsing command-line options\n";
	exit (EXIT_FAILURE);

      } else {
	argvec.push_front (implicitSwitches.front());
	if (implicitSwitches.size() > 1 || !unlimitImplicitSwitches)
	  implicitSwitches.pop_front();
	return true;
      }
    }
  }
  return false;
}

int OptParser::parseUnknownCommand (const string& command, const char* version) {
  if (command == "help" || command == "-help" || command == "--help" || command == "-h") {
    cout << text;
    return EXIT_SUCCESS;
    
  } else if (command == "version" || command == "-version" || command == "--version" || command == "-V") {
    cout << prog << ' ' << version << endl;
    return EXIT_SUCCESS;
    
  } else {
    cerr << briefText << "Unrecognized command: " << command << endl;
  }
  
  return EXIT_FAILURE;
}
