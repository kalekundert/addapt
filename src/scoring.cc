#include <algorithm>
#include <cmath>

#include <boost/algorithm/string.hpp>

extern "C" {
  #include <ViennaRNA/structure_utils.h>
  #include <ViennaRNA/part_func.h>
  #include <ViennaRNA/fold.h>
}

#include <sgrna_design/scoring.hh>

namespace sgrna_design {

ScoreFunction::ScoreFunction() {}

void 
ScoreFunction::operator+=(ScoreTermPtr term) {
	my_terms.push_back(term);
}

double
ScoreFunction::evaluate(ConstructPtr sgrna) const {
	ViennaRnaFold apo_fold(sgrna, LigandEnum::NONE);
	ViennaRnaFold holo_fold(sgrna, LigandEnum::THEO);

	double score = 0;

	for(ScoreTermPtr term: my_terms) {
		score += term->weight() * term->evaluate(sgrna, apo_fold, holo_fold);
	}

	return score;
}


ViennaRnaFold::ViennaRnaFold(ConstructConstPtr sgrna, LigandEnum ligand):

	my_sgrna(sgrna),
	my_seq(sgrna->seq()) {

	// Upper-casing the sequence is critically important!  Without this step, 
	// ViennaRNA will silently produce incorrect results.  I realized I needed to 
	// do this by carefully reading RNAfold.c when I realized that rhf(6) wasn't 
	// being folded correctly.
	boost::to_upper(my_seq);

	my_fc = vrna_fold_compound(
			my_seq.c_str(), NULL, VRNA_OPTION_MFE | VRNA_OPTION_PF);

	my_fold = new char [my_seq.length() + 1];

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

	vrna_pf(my_fc, my_fold);

}

ViennaRnaFold::~ViennaRnaFold() {
	if (my_fc) {
		vrna_fold_compound_free(my_fc);
		delete [] my_fold;
	}
}

char const *
ViennaRnaFold::base_pair_string() const {
	return my_fold;
}

double
ViennaRnaFold::base_pair_prob(int a, int b) const {
	auto indices = normalize_range(my_seq, a, b, IndexEnum::ITEM);
	// The ViennaRNA matrices are 1-indexed.
	int i = indices.first + 1;
	int j = indices.second + 1;
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
		ConditionEnum condition,
		vector<string> selection_a,
		vector<string> selection_b,
		double weight):

	ScoreTerm(weight),
	my_condition(condition),
	my_selection_a(selection_a),
	my_selection_b(selection_b) {}

double
BasePairingTerm::evaluate(
		ConstructConstPtr sgrna,
		RnaFold const &apo_fold,
		RnaFold const &holo_fold) const {

	vector<int> indices_a = domains_to_indices(sgrna, my_selection_a);
	vector<int> indices_b = domains_to_indices(sgrna, my_selection_b);
	RnaFold const & fold =
		(my_condition == ConditionEnum::HOLO) ? holo_fold : apo_fold;

	return prob_paired(indices_a, indices_b, fold);

}


// Score terms I might want
// ========================
// 1. RequiredBasePairing
//
// 2. Rename SpecificLigandSensitivityTerm to ConditionalBasePairing
//
// 3. Conditionally unpaired (i.e. to prevent misfolding the nexus).
//

LigandSensitivityTerm::LigandSensitivityTerm(
		vector<string> selection,
		double weight):

	ScoreTerm(weight),
	my_selection(selection) {

	if(my_selection.empty()) {
		throw "cannot score empty selection";
	}
}

double
LigandSensitivityTerm::evaluate(
		ConstructConstPtr sgrna,
		RnaFold const &apo_fold,
		RnaFold const &holo_fold) const {

	// Return the number of base pairs (combined between the apo and holo states) 
	// weighted by how much more likely each base pair is in one state compared 
	// to the other.

	double sensitivity = 0;
	double p_apo, p_holo;

	vector<int> indices = domains_to_indices(sgrna, my_selection);

	for(auto it_i = indices.begin(); it_i != indices.end(); it_i++) {
		for(auto it_j = it_i; it_j != indices.end(); it_j++) {
			p_apo = apo_fold.base_pair_prob(*it_i, *it_j);
			p_holo = holo_fold.base_pair_prob(*it_i, *it_j);

			if(p_apo > 0 and p_holo > 0) {
				sensitivity += p_apo * log(p_apo / p_holo);
				sensitivity += p_holo * log(p_holo / p_apo);
			}
		}
	}

	// There is some cancellation going on here.  We should double the 
	// normalization factor because we counted base pairs in apo state and the 
	// holo state, but we should halve it because we can only expect N 
	// nucleotides to form at most N/2 base pairs.

	return sensitivity / indices.size();
}


SpecificLigandSensitivityTerm::SpecificLigandSensitivityTerm(
		ConditionEnum condition,
		vector<string> selection,
		vector<string> targets,
		double weight):

	ScoreTerm(weight),
	my_condition(condition),
	my_selection(selection),
	my_targets(targets) {
	
	if(my_selection.empty()) {
		throw "cannot score empty selection";
	}
}

double
SpecificLigandSensitivityTerm::evaluate(
		ConstructConstPtr sgrna,
		RnaFold const &apo_fold,
		RnaFold const &holo_fold) const {

	vector<int> selection_indices = domains_to_indices(sgrna, my_selection);
	vector<int> target_indices = domains_to_indices(sgrna, my_targets);

	double sensitivity = 0;
	double p_apo, p_holo;

	for(int i: selection_indices) {
		for(int j: target_indices) {
			p_apo = apo_fold.base_pair_prob(i, j);
			p_holo = holo_fold.base_pair_prob(i, j);

			if(p_apo > 0 and p_holo > 0) {
				switch(my_condition) {
					case ConditionEnum::APO:
						sensitivity += p_apo * log(p_apo / p_holo);
						break;
					case ConditionEnum::HOLO:
						sensitivity += p_holo * log(p_holo / p_apo);
						break;
				}
			}
		}
	}

	return sensitivity / selection_indices.size();


}

}

namespace std {

ostream &
operator<<(ostream& out, const sgrna_design::ConditionEnum& condition) {
	switch(condition) {
		case sgrna_design::ConditionEnum::APO: out << "APO"; break;
		case sgrna_design::ConditionEnum::HOLO: out << "HOLO"; break;
	}
	return out;
}

}

