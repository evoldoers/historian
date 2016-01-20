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
#define PROPALOR_PROGNAME "propalor"
#define PROPALOR_VERSION  "0.1"

struct PropalorUsage : OptParser {
  PropalorUsage (int argc, char** argv);
};

PropalorUsage::PropalorUsage (int argc, char** argv)
  : OptParser (argc, argv, PROPALOR_PROGNAME, "{help,version} [options]")
{
}

int main (int argc, char** argv) {

  try {
    PropalorUsage usage (argc, argv);
  
    const string command = usage.getCommand();

    return usage.parseUnknownCommand (command, PROPALOR_VERSION);

  } catch (...) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
