#include <sstream>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "logger.h"
#include "regexmacros.h"

Logger logger;

// POSIX basic regular expressions
const regex all_v ("-" RE_PLUS("v"), regex_constants::basic);
const regex numeric_v ("-v" RE_NUMERIC_GROUP, regex_constants::basic);

// functions
string ansiEscape (int code) {
  return string("\x1b[") + to_string(code) + "m";
}

Logger::Logger()
  : verbosity(0), useAnsiColor(true)
{
  for (int col : { 7, 2, 3, 5, 6, 1, 2, 3, 5, 6 })  // no blue, it's invisible
    logAnsiColor.push_back (ansiEscape(30 + col) + ansiEscape(40));
  ansiColorOff = ansiEscape(0);
}

void Logger::addTag (const char* tag) {
  addTag (string (tag));
}

void Logger::addTag (const string& tag) {
  logTags.insert (tag);
}

void Logger::setVerbose (int v) {
  verbosity = max (verbosity, v);
}

void Logger::colorOff() {
  useAnsiColor = false;
}

bool Logger::parseLogArgs (deque<string>& argvec) {
  smatch sm;
  if (argvec.size()) {
    const string& arg = argvec[0];
    if (arg == "-log") {
      Require (argvec.size() > 1, "%s must have an argument", arg.c_str());
      addTag (argvec[1]);
      argvec.pop_front();
      argvec.pop_front();
      return true;

    } else if (arg == "-verbose") {
      setVerbose (1);
      argvec.pop_front();
      return true;

    } else if (regex_match (arg, all_v)) {
      setVerbose ((unsigned int) (arg.size() - 1));
      argvec.pop_front();
      return true;

    } else if (regex_match (arg, sm, numeric_v)) {
      setVerbose (atoi (sm.str(1).c_str()));
      argvec.pop_front();
      return true;

    } else if (arg == "-nocolor") {
      useAnsiColor = false;
      argvec.pop_front();
      return true;
    }
  }
  return false;
}

string Logger::args() const {
  string a;
  if (verbosity > 0)
    a += " -v" + to_string(verbosity);
  for (const auto& t : logTags)
    a += " -log " + t;
  if (!useAnsiColor)
    a += " -nocolor";
  return a;
}


ProgressLogger::ProgressLogger (int verbosity, const char* function, const char* file, int line)
  : msg(NULL), verbosity(verbosity), function(function), file(file), line(line)
{ }

void ProgressLogger::initProgress (const char* desc, ...) {
  startTime = std::chrono::system_clock::now();
  lastElapsedSeconds = 0;
  reportInterval = 2;

  time_t rawtime;
  struct tm * timeinfo;

  time (&rawtime);
  timeinfo = localtime (&rawtime);
  
  va_list argptr;
  va_start (argptr, desc);
  vasprintf (&msg, desc, argptr);
  va_end (argptr);

  if (logger.testVerbosityOrLogTags (verbosity, function, file)) {
    ostringstream l;
    l << msg << ": started at " << asctime(timeinfo);
    logger.print (l.str(), file, line, verbosity);
  }
}

ProgressLogger::~ProgressLogger() {
  if (msg)
    free (msg);
}

void ProgressLogger::logProgress (double completedFraction, const char* desc, ...) {
  va_list argptr;
  const std::chrono::system_clock::time_point currentTime = std::chrono::system_clock::now();
  const auto elapsedSeconds = std::chrono::duration_cast<std::chrono::seconds> (currentTime - startTime).count();
  const double estimatedTotalSeconds = elapsedSeconds / completedFraction;
  if (elapsedSeconds > lastElapsedSeconds + reportInterval) {
    const double estimatedSecondsLeft = estimatedTotalSeconds - elapsedSeconds;
    const double estimatedMinutesLeft = estimatedSecondsLeft / 60;
    const double estimatedHoursLeft = estimatedMinutesLeft / 60;
    const double estimatedDaysLeft = estimatedHoursLeft / 24;

    if (completedFraction > 0 && logger.testVerbosityOrLogTags (verbosity, function, file)) {
      char *progMsg;
      va_start (argptr, desc);
      vasprintf (&progMsg, desc, argptr);
      va_end (argptr);

      ostringstream l;
      l << msg << ": " << progMsg << ". Estimated time left: ";
      if (estimatedDaysLeft > 2)
	l << estimatedDaysLeft << " days";
      else if (estimatedHoursLeft > 2)
	l << estimatedHoursLeft << " hrs";
      else if (estimatedMinutesLeft > 2)
	l << estimatedMinutesLeft << " mins";
      else
	l << estimatedSecondsLeft << " secs";
      l << " (" << (100*completedFraction) << "%)" << endl;

      logger.print (l.str(), file, line, verbosity);
      
      free(progMsg);
    }
    
    lastElapsedSeconds = elapsedSeconds;
    reportInterval = fmin (10., 2*reportInterval);
  }
}
