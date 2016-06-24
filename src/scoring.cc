#include <algorithm>
#include <cmath>

extern "C" {
  #include <ViennaRNA/structure_utils.h>
  #include <ViennaRNA/part_func.h>
  #include <ViennaRNA/fold.h>
}

#include <sgrna_design/scoring.hh>

namespace sgrna_design {

enum class LigandEnum {
	NONE,
	THEO,
};


ScoreFunction::ScoreFunction() {}

void 
ScoreFunction::operator+=(ScoreTermPtr term) {
	my_terms.push_back(term);
}

double
ScoreFunction::evaluate(ConstructPtr sgrna) const {
	RnaFold no_lig(sgrna, LigandEnum::NONE);
	RnaFold lig(sgrna, LigandEnum::THEO);

	double score = 0;

	for(ScoreTermPtr term: my_terms) {
		score += term->weight() * term->evaluate(sgrna, no_lig, lig);
	}

	return score;
}


RnaFold::RnaFold(ConstructPtr sgrna, LigandEnum ligand):

	my_sgrna(sgrna),
	my_seq(sgrna->seq()) {

	my_fc = vrna_fold_compound(
			my_seq.c_str(), NULL, VRNA_OPTION_MFE | VRNA_OPTION_PF);

	char const *aptamer_seq = NULL;
	char const *aptamer_fold = NULL;
	double aptamer_kd_uM = 0;

	switch(ligand) {

		case LigandEnum::NONE :
			break;

		case LigandEnum::THEO :
			aptamer_seq =  "GAUACCAGCCGAAAGGCCCUUGGCAGC";
			aptamer_fold = "(...((.(((....)))....))...)";
			aptamer_kd_uM = 0.32;
			break;
	}

	if (aptamer_seq) {
		//double rt_37 = 1.9858775e-3 * 310;  // kcal/mol at 37°C
		//double std_conc = 1e6;              // 1M in μM
		double kT = my_fc->exp_params->kT / 1000;
		double aptamer_dg = kT * log(aptamer_kd_uM / 1e6);
		vrna_sc_add_hi_motif(
				my_fc, aptamer_seq, aptamer_fold, aptamer_dg, VRNA_OPTION_DEFAULT);
	}

	vrna_pf(my_fc, NULL);

}

RnaFold::~RnaFold() {
	if (my_fc) {
		vrna_fold_compound_free(my_fc);
	}
}

double
RnaFold::base_pair_prob(Nucleotide nuc_a, Nucleotide nuc_b) const {
	return base_pair_prob(my_sgrna->index(nuc_a), my_sgrna->index(nuc_b));
}

double
RnaFold::base_pair_prob(Nucleotide nuc_a, int b) const {
	return base_pair_prob(my_sgrna->index(nuc_a), b);
}

double
RnaFold::base_pair_prob(int a, int b) const {
	int i = std::min(a, b) + 1;
	int j = std::max(a, b) + 1;
	return (i != j) ? my_fc->exp_matrices->probs[my_fc->iindx[i] - j] : 0;
}
	

ScoreTerm::ScoreTerm() : ScoreTerm(0) {}

ScoreTerm::ScoreTerm(double weight) : my_weight(weight) {}

double
ScoreTerm::weight() const {
	return my_weight;
}

void
ScoreTerm::weight(double weight) {
	my_weight = weight;
}


double
prob_paired(
		vector<int> indices_i,
		vector<int> indices_j,
		RnaFold const & fold) {

	double prob = 0;

	for(int i: indices_i) {
		for(int j: indices_j) {
			prob += fold.base_pair_prob(i, j);
		}
	}

	return prob;
}

vector<int>
domains_to_indices(
		ConstructConstPtr sgrna,
		vector<string> domains) {

	vector<int> indices;

	for(string domain: domains) {
		for(int i = sgrna->index_5(domain); i <= sgrna->index_3(domain); i++) {
			indices.push_back(i);
		}
	}

	return indices;
}


BasePairingTerm::BasePairingTerm(
		vector<string> selection,
		vector<string> no_lig_targets,
		vector<string> lig_targets,
		double weight):

	ScoreTerm(weight),
	my_selection(selection),
	my_no_lig_targets(no_lig_targets),
	my_lig_targets(lig_targets) {}

double
BasePairingTerm::evaluate(
		ConstructConstPtr sgrna,
		RnaFold const &no_lig_fold,
		RnaFold const &lig_fold) const {

	vector<int> sele_indices = domains_to_indices(sgrna, my_selection);
	vector<int> no_lig_indices = domains_to_indices(sgrna, my_no_lig_targets);
	vector<int> lig_indices = domains_to_indices(sgrna, my_lig_targets);

	double no_lig_base_pairing = prob_paired(
			sele_indices, no_lig_indices, no_lig_fold);
	double lig_base_pairing = prob_paired(
			sele_indices, lig_indices, lig_fold);

	return no_lig_base_pairing * lig_base_pairing;

}


}

