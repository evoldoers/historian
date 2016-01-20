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
#define IDHIST_PROGNAME "idhist"
#define IDHIST_VERSION  "0.1"

struct ProgUsage : OptParser {
  ProgUsage (int argc, char** argv);
};

ProgUsage::ProgUsage (int argc, char** argv)
  : OptParser (argc, argv, IDHIST_PROGNAME, "{help,version} [options]")
{
}

int main (int argc, char** argv) {

  try {
    ProgUsage usage (argc, argv);
  
    const string command = usage.getCommand();

    return usage.parseUnknownCommand (command, IDHIST_VERSION);

  } catch (...) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
