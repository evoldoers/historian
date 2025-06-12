#ifndef LOGGER_INCLUDED
#define LOGGER_INCLUDED

#include <list>
#include <set>
#include <map>
#include <string>
#include <deque>
#include <ratio>
#include <chrono>
#include <iostream>
#include <sstream>
#include "util.h"
#include "vguard.h"

using namespace std;

class Logger {
private:
  int verbosity;
  set<string> logTags;
  bool useAnsiColor;
  vguard<string> logAnsiColor;
  string ansiColorOff;
  
public:
  Logger();
  // configuration
  void addTag (const char* tag);
  void addTag (const string& tag);
  void setVerbose (int v);
  void colorOff();
  bool parseLogArgs (deque<string>& argvec);
  string args() const;
  
  inline bool testVerbosity (int v) {
    return verbosity >= v;
  }

  inline bool testLogTag (const char* tag) {
    return logTags.find(tag) != logTags.end();
  }

  inline bool testVerbosityOrLogTags (int v, const char* tag1, const char* tag2) {
    return verbosity >= v || testLogTag(tag1) || testLogTag(tag2);
  }

  template<class T>
  void print (const T& t, const char* file, int line, int v) {
    if (useAnsiColor) {
      if (v >= 0 && v < logAnsiColor.size())
        clog << logAnsiColor[v];
      else
        clog << logAnsiColor[logAnsiColor.size() - 1];
    }
    clog << t;
    if (useAnsiColor)
      clog << ansiColorOff;
  }
};

extern Logger logger;

#define LoggingAt(V)     (logger.testVerbosity(V))
#define LoggingThisAt(V) (logger.testVerbosityOrLogTags(V,__func__,__FILE__))
#define LoggingTag(T)    (logger.testLogTag(T))

#define LogStream(V,S) do { ostringstream tmpLog; tmpLog << S; logger.print(tmpLog.str(),__FILE__,__LINE__,V); } while(0)

#define LogAt(V,S)     do { if (LoggingAt(V)) LogStream(V,S); } while(0)
#define LogThisAt(V,S) do { if (LoggingThisAt(V)) LogStream(V,S); } while(0)
#define LogThisIf(X,S) do { if (X) LogStream(0,S); } while(0)


/* progress logging */
class ProgressLogger {
public:
  std::chrono::system_clock::time_point startTime;
  double lastElapsedSeconds, reportInterval;
  char* msg;
  int verbosity;
  const char *function, *file;
  int line;
  ProgressLogger (int verbosity, const char* function, const char* file, int line);
  ~ProgressLogger();
  void initProgress (const char* desc, ...);
  void logProgress (double completedFraction, const char* desc, ...);
private:
  ProgressLogger (const ProgressLogger&) = delete;
  ProgressLogger& operator= (const ProgressLogger&) = delete;
};

#define ProgressLog(PLOG,V) ProgressLogger PLOG (V, __func__, __FILE__, __LINE__)

#endif /* LOGGER_INCLUDED */

