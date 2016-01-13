#ifndef UTIL_INCLUDED
#define UTIL_INCLUDED

#include <numeric>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <functional>
#include <cassert>
#include <mutex>
#include <sys/stat.h>

/* uncomment to enable NaN checks */
#define NAN_DEBUG

/* Errors, warnings, assertions.
   Fail(...) and Require(...) are quieter versions of Abort(...) and Assert(...)
   that do not print a stack trace or throw an exception,
   but merely call exit().
   Test(...) does not exit or throw an exception,
   just prints a warning and returns false if the assertion fails.
   Desire(...) is a macro wrapper for Test(...)
   that returns false from the calling function if the test fails.
*/
void Abort(const char* error, ...);
void Assert(int assertion, const char* error, ...);
void Warn(const char* warning, ...);
void Fail(const char* error, ...);
void Require(int assertion, const char* error, ...);
bool Test(int assertion, const char* error, ...);
#define Desire(...) do { if (!Test(__VA_ARGS__)) return false; } while (0)

/* singular or plural? */
std::string plural (long n, const char* singular);
std::string plural (long n, const char* singular, const char* plural);

/* stringify */
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

/* pipe to string */
std::string pipeToString (const char* command, int* status = NULL);

/* temp files and dirs */
#define DefaultTempFilePrefix "tempfile"
#define DefaultTempDirPrefix  "tempdir"

/* temp file */
class TempFile {
private:
  static std::string makeNewPath (std::string basePath, bool usePid);
public:
  static const std::string dir;  /* /tmp */
  static unsigned int count;
  static std::mutex mx;
  static std::string newPathWithPid (std::string basePath);
  static std::string newPath (std::string basePath);
  std::string fullPath;
  TempFile();
  TempFile (const std::string& contents, const char* filenamePrefix = DefaultTempFilePrefix);
  ~TempFile();
  void init (const std::string& contents, const char* filenamePrefix = DefaultTempFilePrefix);
};

/* temp directory */
struct TempDir {
  std::string fullPath;
  bool cleanup;
  TempDir (const char* filenamePrefix = DefaultTempDirPrefix);  // by default, creates tempdir in cwd
  ~TempDir();
  void init();
  std::string makePath (const std::string& filename) const;
};

/* test file exists */
bool file_exists (const char *filename);

/* join */
template<class Container>
std::string join (const Container& c, const char* sep = " ") {
  std::string j;
  for (const auto& s : c) {
    if (!j.empty())
      j += sep;
    j += s;
  }
  return j;
}

/* join */
template<class Container>
std::string to_string_join (const Container& c, const char* sep = " ") {
  std::ostringstream j;
  int n = 0;
  for (const auto& s : c) {
    if (n++ > 0)
      j << sep;
    j << s;
  }
  return j.str();
}

/* sgn function
   http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
 */
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

/* index sort
   http://stackoverflow.com/questions/10580982/c-sort-keeping-track-of-indices
 */
template <typename T>
std::vector<size_t> orderedIndices (std::vector<T> const& values) {
    std::vector<size_t> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<size_t>(0));

    std::sort(
        begin(indices), end(indices),
        [&](size_t a, size_t b) { return values[a] < values[b]; }
    );
    return indices;
}

/* vector sum */
template <typename T>
std::vector<T> vector_sum(const std::vector<T>& a, const std::vector<T>& b)
{
  assert(a.size() == b.size());

  std::vector<T> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), b.begin(), 
		 std::back_inserter(result), std::plus<T>());
  return result;
}

/* vector-scalar product */
template <typename T>
std::vector<T> vector_scale(const T x, const std::vector<T>& a)
{
  std::vector<T> result = a;
  for (auto& y : result)
    y *= x;

  return result;
}

/* escaping a string
   http://stackoverflow.com/questions/2417588/escaping-a-c-string
 */
template<class OutIter>
OutIter write_quoted_escaped(std::string const& s, OutIter out) {
  *out++ = '"';
  for (std::string::const_iterator i = s.begin(), end = s.end(); i != end; ++i) {
    unsigned char c = *i;
    if (' ' <= c and c <= '~' and c != '\\' and c != '"') {
      *out++ = c;
    }
    else {
      *out++ = '\\';
      switch(c) {
      case '"':  *out++ = '"';  break;
      case '\\': *out++ = '\\'; break;
      case '\t': *out++ = 't';  break;
      case '\r': *out++ = 'r';  break;
      case '\n': *out++ = 'n';  break;
      default:
        char const* const hexdig = "0123456789ABCDEF";
        *out++ = 'x';
        *out++ = hexdig[c >> 4];
        *out++ = hexdig[c & 0xF];
      }
    }
  }
  *out++ = '"';
  return out;
}
    
#endif /* UTIL_INCLUDED */
