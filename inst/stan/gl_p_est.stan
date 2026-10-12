functions{
  vector low_coverage_whole_genome_obs(vector empty_vector_paramater, vector p_slice, array[] real BetaBin_theta_slice, array[] int Group_slice){
    int N_patch = size(p_slice), G2;
    real lq = 0;
    int pos = 1;
    vector [2*max(Group_slice)+1] linspace = linspaced_vector(2*max(Group_slice)+1, 0.0, 2.0*max(Group_slice));

    for (i in 1:N_patch){
      // reverse(linspace)*log(1-p) + linspace*log(p) = 2G*log(1-p) + linspace*logit(p)
      // Factor 2G*log(1-p) out of log_sum_exp as a scalar.
      if (Group_slice[i] > 0) {
        G2 = 2*Group_slice[i];
        lq += G2 * log1m(p_slice[i]) + log_sum_exp(to_vector(BetaBin_theta_slice[pos:(pos+G2)]) + linspace[1:(G2+1)] * logit(p_slice[i]));
      }
      pos += 2*Group_slice[i]+1;
    }
    return [lq]';
  }
}

data {
  int < lower = 1 > N_patch ; // number of patches
  int < lower = 1 > N_locus ; // number of loci
  array [N_patch] int N_table; // individuals in each patch
  array [N_locus, 2*sum(N_table) + N_patch] real BetaBin_theta; // genotype liklihoods
  array [N_locus, N_patch] int <lower = 0>  Group; // number of individulas with genotype liklihoods in each patch. Always less than equal to N_table
}

parameters {
  array [N_locus] vector <lower = 0, upper = 1> [N_patch] p;
  vector [0] empty_vector_paramater;
}

model {
  target += sum(map_rect(low_coverage_whole_genome_obs, empty_vector_paramater, p, BetaBin_theta, Group));
}
