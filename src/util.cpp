// feature test macro requirement for ftw
// #define _XOPEN_SOURCE 500

// includes
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <ftw.h>

#include "util.h"
#include "stacktrace.h"
#include "logger.h"

// buffer size for popen
#define PIPE_BUF_SIZE 1024

// buffer size for getcwd
#define DIR_BUF_SIZE 4096


// recursive rmdir
// http://stackoverflow.com/questions/3184445/how-to-clear-directory-contents-in-c-on-linux-basically-i-want-to-do-rm-rf/3184915#3184915
int rmdirCallback(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf)
{
    int rv = remove(fpath);

    if (rv)
        perror(fpath);

    return rv;
}

int rmdirRecursive(const char *path)
{
    return nftw(path, rmdirCallback, 64, FTW_DEPTH | FTW_PHYS);
}

// function defs
void Warn(const char* warning, ...) {
  va_list argptr;
  fprintf(stderr,"Warning: ");
  va_start (argptr, warning);
  vfprintf(stderr,warning,argptr);
  fprintf(stderr,"\n");
  va_end (argptr);
}

void Abort(const char* error, ...) {
  va_list argptr;
  va_start (argptr, error);
  fprintf(stderr,"Abort: ");
  vfprintf(stderr,error,argptr);
  fprintf(stderr,"\n");
  va_end (argptr);
  printStackTrace();
  throw;
}

void Assert(int assertion, const char* error, ...) {
  va_list argptr;
  if(!assertion) {
    va_start (argptr, error);
    fprintf(stderr,"Assertion Failed: ");
    vfprintf(stderr,error,argptr);
    fprintf(stderr,"\n");
    va_end (argptr);
    printStackTrace();
    throw;
  }
}

void Fail(const char* error, ...) {
  va_list argptr;
  va_start (argptr, error);
  vfprintf(stderr,error,argptr);
  fprintf(stderr,"\n");
  va_end (argptr);
  exit (EXIT_FAILURE);
}

void Require(int assertion, const char* error, ...) {
  va_list argptr;
  if(!assertion) {
    va_start (argptr, error);
    vfprintf(stderr,error,argptr);
    fprintf(stderr,"\n");
    va_end (argptr);
    exit (EXIT_FAILURE);
  }
}

bool Test(int assertion, const char* error, ...) {
  va_list argptr;
  if(!assertion) {
    va_start (argptr, error);
    vfprintf(stderr,error,argptr);
    fprintf(stderr,"\n");
    va_end (argptr);
  }
  return assertion;
}

std::string plural (long n, const char* singular) {
  std::string s = std::to_string(n) + " " + singular;
  if (n != 1)
    s += "s";
  return s;
}

std::string plural (long n, const char* singular, const char* plural) {
  std::string s = std::to_string(n) + " " + (n == 1 ? singular : plural);
  return s;
}

const string TempFile::dir = "/tmp";
unsigned int TempFile::count = 0;
std::mutex TempFile::mx;

std::string TempFile::makeNewPath (std::string basePath, bool usePid) {
  std::string fullPath;
  mx.lock();
  do {
    fullPath = basePath + (usePid ? (std::to_string(getpid()) + '.') : std::string()) + std::to_string(++count);
  } while (file_exists (fullPath.c_str()));
  mx.unlock();
  return fullPath;
}

std::string TempFile::newPathWithPid (std::string basePath) {
  return makeNewPath (basePath, true);
}

std::string TempFile::newPath (std::string basePath) {
  return makeNewPath (basePath, false);
}

TempFile::TempFile (const std::string& contents, const char* filenamePrefix) {
  init (contents, filenamePrefix);
}

void TempFile::init (const std::string& contents, const char* filenamePrefix) {
  fullPath = newPathWithPid (dir + filenamePrefix);
  std::ofstream out (fullPath);
  Assert (out.is_open() && !out.fail(), "Couldn't write to temp file %s", fullPath.c_str());
  out << contents;
}

TempFile::~TempFile() {
  if (fullPath.size())
    unlink (fullPath.c_str());
}

TempDir::TempDir (const char* filenamePrefix)
  : cleanup(false)
{
  char dirbuf[DIR_BUF_SIZE];
  fullPath = TempFile::newPath (string (getcwd (dirbuf, DIR_BUF_SIZE)) + '/' + filenamePrefix);
}

void TempDir::init() {
  if (!file_exists(fullPath.c_str())) {
    Assert (mkdir(fullPath.c_str(),0700) == 0, "Couldn't make temp directory %s", fullPath.c_str());
    cleanup = true;
  }
}

TempDir::~TempDir() {
  if (cleanup && fullPath.size() && file_exists(fullPath.c_str()))
    rmdirRecursive (fullPath.c_str());
}

string TempDir::makePath (const string& filename) const {
  return fullPath + '/' + filename;
}

string pipeToString (const char* command, int* status) {
  string result;
  FILE* pipe = popen (command, "r");
  char line[PIPE_BUF_SIZE];

  while (fgets(line, PIPE_BUF_SIZE, pipe))
    result += line;

  const int s = pclose (pipe);
  if (status)
    *status = s;
  
  return result;
}

bool file_exists (const char *filename) {
  struct stat buffer;   
  return stat (filename, &buffer) == 0;
}
