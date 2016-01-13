#include "../src/kseq.h"
#include "../src/fastseq.h"
#include "../src/gason.h"
#include "../src/jsonutil.h"
#include "../src/knhx.h"
#include "../src/logger.h"
#include "../src/logsumexp.h"
#include "../src/vguard.h"
#include "../src/optparser.h"

// GNU --version
#define PROPAELOR_PROGNAME "propaelor"
#define PROPAELOR_VERSION  "0.1"

struct PropaelorUsage : OptParser {
  PropaelorUsage (int argc, char** argv);
};

PropaelorUsage::PropaelorUsage (int argc, char** argv)
  : OptParser (argc, argv, PROPAELOR_PROGNAME, "{help,version} [options]")
{
}

int main (int argc, char** argv) {

  try {
    PropaelorUsage usage (argc, argv);
  
    const string command = usage.getCommand();

    return usage.parseUnknownCommand (command, PROPAELOR_VERSION);

  } catch (...) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
