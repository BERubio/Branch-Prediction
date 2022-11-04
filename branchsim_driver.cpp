#include <cstdio>
#include <cinttypes>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include "branchsim.hpp"

void print_help_and_exit(void) {
  printf("branchsim [OPTIONS] < traces/file.trace\n");
  printf("  -k [NUM_ENTRIES]\tNumber of entries in the predictor table\n");
  printf("  -c [CBITS]\t\tNumber of bits for each counter\n");
  printf("  -s [HBITS]\t\tNumber of bits for the history register (if needed)\n");
  printf("  -p [PTYPE]\t\tType of predictor to use (S, L, or G)\n");
  printf("  -h\t\t\tThis helpful output\n");

  exit(0);
}

void print_statistics(branch_stats_t* p_stats);

int main(int argc, char* argv[]) {
  int opt;
  predictor_type ptype= PTYPE_BIMODAL;
  int num_entries     = 1;
  int counter_bits    = 1;
  int history_bits    = 1;

  // Process arguments
  while(-1 != (opt = getopt(argc, argv, "k:c:s:p:h"))) {
    switch(opt) {
      case 'k':
        num_entries = atoi(optarg);
        break;
      case 'c':
        counter_bits = atoi(optarg);
        break;
      case 's':
        history_bits = atoi(optarg);
        break;
      case 'p':
        ptype = static_cast<predictor_type>(optarg[0]);
        break;
      case 'h':
        // Fall through
      default:
        print_help_and_exit();
        break;
    }
  }

  printf("Branch Predictor Settings\n");
  printf("Predictor Type: %c\n", static_cast<char>(ptype));
  printf("# Entries: %d\n", num_entries);
  printf("# Bits Per Counter: %d\n", counter_bits);
  printf("# Bits Per History: %d\n", history_bits);
  printf("\n");

  branch_stats_t stats;
  memset(&stats, 0, sizeof(branch_stats_t));

  // Do some setup
  setup_predictor(ptype, num_entries, counter_bits, history_bits, &stats);

  // For each trace in the file
  while (!feof(stdin)) {
    std::uint64_t pc = 0;
    char result = 0;
    int ret = fscanf(stdin, "%" PRIx64 " %c\n", &pc, &result);

    if(ret == 2) {
      branch_dir actual = static_cast<branch_dir>(result);
      branch_dir predicted = predict_branch(pc, &stats);

      update_predictor(pc, actual, predicted, &stats);
    }
  }
  complete_predictor(&stats);

  print_statistics(&stats);

  return 0;
}

void print_statistics(branch_stats_t* p_stats) {
  printf("Branch Predictor Statistics\n");
  printf("# Branches: %" PRIu64 "\n", p_stats->num_branches);
  printf("# Predicted Taken: %" PRIu64 "\n",   p_stats->pred_taken);
  printf("# Predicted Not Taken: %" PRIu64 "\n",   p_stats->pred_not_taken);
  printf("# Correct: %" PRIu64 "\n",   p_stats->correct);
  printf("Misprediction Rate: %f\n",   p_stats->misprediction_rate);
  printf("Storage overhead: %" PRIu64 " bits\n",   p_stats->storage_overhead);
}

