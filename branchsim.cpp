#include "branchsim.hpp"

/**
 * A set of global variables
 * You may define and set any global classes and variables here as needed.
 */
static unsigned int* pht;
static std::uint64_t ghr;
static unsigned int ctr_max;
static unsigned int ctr_weakly_taken;
static std::uint64_t index_mask;
static std::uint64_t  history_shift;
static std::uint64_t history_mask;
static predictor_type pred_type;
static std::uint64_t* hist_table;

/**
 * This subroutine initializes the branch predictors.
 * You many add and initialize any global or heap
 * variables as needed.
 * (1) You're responsible for completing this routine
 *
 * Inputs and outputs:
 * @param[in]   ptype       The type of branch predictor to simulate
 * @param[in]   num_entries The number of entries a PC is hashed into
 * @param[in]   counter_bits The number of bits per counter
 * @param[in]   history_bits The number of bits per history
 * @param[out]  p_stats     Pointer to the stats structure
 */
void setup_predictor(predictor_type ptype,
                     int num_entries,
                     int counter_bits,
                     int history_bits,
                     branch_stats_t* p_stats) {
  ctr_max = (unsigned int)(1 << counter_bits) - 1;
  ctr_weakly_taken = (unsigned int) (1 << (counter_bits - 1));

  index_mask = (std::uint64_t)num_entries - 1;
  history_mask = (std::uint64_t)(1 << history_bits) - 1;
  history_shift = history_bits;

  int pht_size = num_entries;
  pred_type = ptype;
  /*
   * For each type of branch prediction, complete the following switch
   * to define pht_size and p_stats based on num_entries, counter_bits, and history_bits.
   * The code for Bimodal branch predictor is provided as an example.
   */
  switch (ptype) {
    case PTYPE_BIMODAL: {
      //num_entries = 2^n 
      //counter_bits = m-bit counter per entry in pht
      pht_size = num_entries;
      p_stats->storage_overhead = ((std::uint64_t)pht_size) * ((std::uint64_t)counter_bits);
      break;
    }
    case PTYPE_TWO_LEVEL_ADAPTIVE: {
      //pht_size is determined by the num_entries in the HISTORY table NOT PHT
      pht_size = (1 << (history_bits+1));
      p_stats->storage_overhead = (((std::uint64_t)history_bits) * ((std::uint64_t)num_entries)) + (((std::uint64_t)pht_size));
      break;
    }
    case PTYPE_LOCAL_HISTORY: {
      //pht_size is determined by the num_entries in the HISTORY table NOT PHT
      pht_size = (1 << history_bits); 
      p_stats->storage_overhead = (((std::uint64_t)history_bits) * ((std::uint64_t)num_entries)) + ((((std::uint64_t)counter_bits)*((std::uint64_t)pht_size)));
      break;
    }
    case PTYPE_GSHARE: {
      pht_size = num_entries * counter_bits;
      //How do I incorporate GHR (n-bits) into storage overhead???
      p_stats->storage_overhead = (((std::uint64_t)pht_size) + ((std::uint64_t)history_bits));
      break;
    }
  }

  /*
   * We initial the pht, ghr, and hist_table in the following:
   */
  pht = new unsigned int[pht_size];
  for (int i = 0; i < pht_size; i++) {
    pht[i] = ctr_weakly_taken;
  }
  ghr = 0;
  if (pred_type == PTYPE_LOCAL_HISTORY || pred_type == PTYPE_TWO_LEVEL_ADAPTIVE) {
    hist_table = new std::uint64_t[num_entries];
    for (int i = 0; i < num_entries; i++) {
      hist_table[i] = 0;
    }
  }
}

/**
 * This subroutine calculates the index for accessing the pht
 * (2) You're responsible for completing this routine
 *
 * Inputs and outputs:
 * @param[in]   pc       The program counter
 * @param[out]  The index for accessing the pht
 * The code for bimodal branch predictor is provided as an example.
 */
std::uint64_t get_index(std::uint64_t pc) {
  /*
   * For each type of branch prediction, complete the following switch
   * to calculate the index. You may combine pc with index_mask and use
   * history_shift, hist_table, and history_mask in your calculation.
   * The code for the Local History branch predictor is provided as an example.
   */
  switch (pred_type) {
    case PTYPE_BIMODAL: {
      std::uint64_t hist_index = pc;
      return (hist_index * (1 << history_shift));
    }
    case PTYPE_TWO_LEVEL_ADAPTIVE: {
      //Calculating the index should be similar to LOCAL minus the fact that there is only one PHT
      std::uint64_t hist_index = pc & index_mask;
      return (hist_table[hist_index] & history_mask);
    }
    case PTYPE_LOCAL_HISTORY: {
      //LOCAL_HISTORY get_index code provided as an example
      std::uint64_t hist_index = pc & index_mask;
      return (hist_index * (1 << history_shift)) + (hist_table[hist_index] & history_mask);
    }
    case PTYPE_GSHARE: {
      //This code is generating a seg. fault...FIXED: Had to mask ghr with index mask before XOR-ing
      std::uint64_t hist_index = pc & index_mask;
      return hist_index ^ (ghr & index_mask);
    }
    default:
      return 0;
  }
}

/**
 * This subroutine run the branch prediction for a PC, returns either TAKEN or
 * NOT_TAKEN and accordingly increaments the pred_taken or pred_not_taken++.
 * (3) You're responsible for completing this routine
 *
 * @param[in]   pc          The PC value of the branch instruction.
 * @param[out]  p_stats     Pointer to the stats structure
 *
 * @return                  Either TAKEN ('T'), or NOT_TAKEN ('N')
 */
branch_dir predict_branch(std::uint64_t pc, branch_stats_t* p_stats) {
  // Increment branch count
  p_stats->num_branches++;

  // Identify index
  std::uint64_t index = get_index(pc);

  // Predict the branch by accessing the pht
  switch (pred_type) {
    // predcit for bimodal: get_index(pc)= index -> pht[index] -> prediction
    case PTYPE_BIMODAL: {
      if(pht[index] >= ctr_weakly_taken){
        p_stats->pred_taken++;
        return TAKEN;
      }else{
        p_stats->pred_not_taken++;
        return NOT_TAKEN;
      }
      break;
    }
    case PTYPE_TWO_LEVEL_ADAPTIVE: {
      if(pht[index] >= ctr_weakly_taken){
        p_stats->pred_taken++;
        return TAKEN;
      }else{
        p_stats->pred_not_taken++;
        return NOT_TAKEN;
      }
      break;
    }
    //Local history prediction feels different than others...May need different implementation?
    case PTYPE_LOCAL_HISTORY: {
      if(pht[index]>=ctr_weakly_taken){
        p_stats->pred_taken++;
        return TAKEN;
      }else{
        p_stats->pred_not_taken++;
        return NOT_TAKEN;
      }
      break;
    }
    case PTYPE_GSHARE: {
      if(pht[index] >= ctr_weakly_taken){
        p_stats->pred_taken++;
        return TAKEN;
      }else{
        p_stats->pred_not_taken++;
        return NOT_TAKEN;
      }
    }
    default:
      return TAKEN;
  }
}

/**
 * This subroutine updates the branch predictor.
 * The branch predictor needs to be notified whether
 * the prediction was right or wrong.
 * To do so, update pht, ghr, hist_table, and any other variables as needed.
 * (4) You're responsible for completing this routine
 *
 * @param[in]   pc          The PC value of the branch instruction.
 * @param[in]   actual      The actual direction of the branch
 * @param[in]   predicted   The predicted direction of the branch (from predict_branch)
 * @param[out]  p_stats     Pointer to the stats structure
 */
void update_predictor(std::uint64_t pc,
                      branch_dir actual,
                      branch_dir predicted,
                      branch_stats_t* p_stats) {
  std::uint64_t index = get_index(pc);

  
  switch (pred_type) {
    //For Bimodal update, just update pht at the index and p_stats
    case PTYPE_BIMODAL: {
      //IF CORRECT:------------
      if (actual == predicted) {
        p_stats->correct++;
        if(actual == TAKEN){
          pht[index] = pht[index]<<1;
        }else{
          pht[index] = pht[index]<<0;
        }
      }else{
        //IF INCORRECT:-----------
        if(actual == TAKEN){
          p_stats->misprediction_rate++;
          pht[index] = pht[index]<<1;
        }else{
          p_stats->misprediction_rate++;
          pht[index] = pht[index]<<0;
        }
      }
      break;
    }
    //For Two-Level, update PHT, history table and P_stats
    case PTYPE_TWO_LEVEL_ADAPTIVE: {
      //IF CORRECT: ---------------
      if(actual == predicted){
        p_stats->correct++;
        if(actual == TAKEN){
          pht[index] = pht[index]+1;
          hist_table[index] = hist_table[index] << 1;
        }else{
          pht[index] = pht[index]-1;
          hist_table[index] = hist_table[index] << 0;
        }
      }else{
        //IF INCORRECT: ------------------
        if(actual== TAKEN){
          p_stats->misprediction_rate++;
          pht[index] = pht[index]+1;
          hist_table[index] = hist_table[index] << 1;
        }else{
          p_stats->misprediction_rate++;
          pht[index] = pht[index]-1;
          hist_table[index] = hist_table[index] << 0;
        }
      }
      break;
    }
    //For Local, update PHT, history table, and P-stats
    case PTYPE_LOCAL_HISTORY: {
      //IF CORRECT:----------
      if(actual == predicted){
        p_stats->correct++;
        if(actual == TAKEN){
          pht[index] = pht[index]+1;
          hist_table[index] = hist_table[index] << 1;
        }else{
          pht[index] = pht[index]-1;
          hist_table[index] = hist_table[index] << 0;
        }
        break;
      }else{
        //IF INCORRECT:--------------
        if(actual== TAKEN){
          p_stats->misprediction_rate++;
          pht[index] = pht[index]+1;
          hist_table[index] = hist_table[index] << 1;
        }else{
          p_stats->misprediction_rate++;
          pht[index] = pht[index]-1;
          hist_table[index] = hist_table[index] << 0;
        }
      }
    }
    case PTYPE_GSHARE: {
      if(actual == predicted){
        if (actual == TAKEN) {
          ghr = (ghr << 1) | 1;
          p_stats->misprediction_rate++;
          pht[index] = pht[index] + 1;
        } else {
          ghr = (ghr << 1) | 0;
          p_stats->misprediction_rate++;
          pht[index] = pht[index] - 1;
        }
      }else{
        if (actual == TAKEN) {
          ghr = (ghr << 1) | 1;
          p_stats->misprediction_rate++;
          pht[index] = pht[index] + 1;
        } else {
          ghr = (ghr << 1) | 0;
          p_stats->misprediction_rate++;
          pht[index] = pht[index] - 1;
        }
      }
      break;
    }
  }
}

/**
 * This subroutine cleans up any outstanding memory operations and calculating overall statistics.
 *
 * @param[out]  p_stats     Pointer to the statistics structure
 */
void complete_predictor(branch_stats_t *p_stats) {
  p_stats->misprediction_rate = (double)(p_stats->num_branches - p_stats->correct) / (double)p_stats->num_branches;
  delete[] pht;
  if (pred_type == PTYPE_LOCAL_HISTORY) {
    delete[] hist_table;
  }
}
