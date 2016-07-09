#include <algorithm>
#include <cmath>

#include <boost/algorithm/string.hpp>

extern "C" {
  #include <ViennaRNA/structure_utils.h>
  #include <ViennaRNA/part_func.h>
  #include <ViennaRNA/fold.h>
	#include <ViennaRNA/constraints.h>
}

#include <sgrna_design/scoring.hh>

namespace sgrna_design {

ViennaRnaFold::ViennaRnaFold(ConstructConstPtr sgrna, LigandEnum ligand):

	my_sgrna(sgrna),
	my_seq(sgrna->seq()),
	my_ligand(ligand),
	my_bppm_fc(nullptr) {

	// Upper-casing the sequence is critically important!  Without this step, 
	// ViennaRNA will silently produce incorrect results.  I realized I needed to 
	// do this by carefully reading RNAfold.c when I realized that rhf(6) wasn't 
	// being folded correctly.
	boost::to_upper(my_seq);
}

ViennaRnaFold::~ViennaRnaFold() {
	for(auto fc: my_fcs) {
		vrna_fold_compound_free(fc);
	}
}

double
ViennaRnaFold::base_pair_prob(int a, int b) const {
	// Perform the partition function calculation if this is the first time a 
	// base-pair probability is being requested.  Cache the result.
	if(my_bppm_fc == nullptr) {
		my_bppm_fc = make_fold_compound(true);
		vrna_pf(my_bppm_fc, NULL);
	}

	auto indices = normalize_range(my_seq, a, b, IndexEnum::ITEM);
	int i = indices.first + 1; // The ViennaRNA matrices are 1-indexed.
	int j = indices.second + 1;
	double *bppm = my_bppm_fc->exp_matrices->probs;
	return (i != j) ? bppm[my_bppm_fc->iindx[i] - j] : 0;
}
	
double
ViennaRnaFold::macrostate_prob(string macrostate) const {
	vrna_fold_compound_t *fc = make_fold_compound(false);

	// Calculate the free energy for the whole ensemble.
	double g_tot = vrna_pf(fc, NULL);

	// Add a constraint that defines the "active" macrostate.
	vrna_constraints_add(fc, macrostate.c_str(), VRNA_CONSTRAINT_DB_DEFAULT);

	// Calculate the free energy for the "active" macrostate.
	double g_active = vrna_pf(fc, NULL);

	// Return the probability that the construct will be in the "active" 
	// macrostate at equilibrium.
	double kT = fc->exp_params->kT / 1000;
	return exp((g_tot - g_active) / kT);
}

vrna_fold_compound_t *
ViennaRnaFold::make_fold_compound(bool compute_bppm) const {
	// Tell ViennaRNA not to calculate the base-pair probability matrix (BPPM) if 
	// we won't be using it.
	vrna_md_t md;
	vrna_md_set_default(&md);
  md.backtrack = compute_bppm;
	md.compute_bpp = compute_bppm;

	// Create a new "fold compound" data structure and store a pointer to it so 
	// it can be deallocated later.
	vrna_fold_compound_t *fc =
		vrna_fold_compound(my_seq.c_str(), &md, VRNA_OPTION_PF);
	my_fcs.push_back(fc);

	// Add the aptamer motif.
	char const *aptamer_seq = NULL;
	char const *aptamer_fold = NULL;
	double aptamer_kd_uM = 0;

	switch(my_ligand) {
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
		double kT = fc->exp_params->kT / 1000;
		double aptamer_dg = kT * log(aptamer_kd_uM / 1e6);
		vrna_sc_add_hi_motif(
				fc, aptamer_seq, aptamer_fold, aptamer_dg, VRNA_OPTION_DEFAULT);
	}

	return fc;
}


ScoreFunction::ScoreFunction() {}

double
ScoreFunction::evaluate(ConstructConstPtr sgrna) const {
	EvaluatedScoreFunction table;
	return evaluate(sgrna, table);
}
		
double
ScoreFunction::evaluate(
		ConstructConstPtr sgrna,
		EvaluatedScoreFunction &table) const {

	ViennaRnaFold apo_fold(sgrna, LigandEnum::NONE);
	ViennaRnaFold holo_fold(sgrna, LigandEnum::THEO);

	double score = 0;
	double evaluated_term;

	table.clear();

	for(ScoreTermPtr term: my_terms) {
		evaluated_term = term->evaluate(sgrna, apo_fold, holo_fold);
		score += term->weight() * evaluated_term;
		table.push_back({term->name(), term->weight(), evaluated_term});
	}

	return score;
}

void 
ScoreFunction::add_term(ScoreTermPtr term) {
	my_terms.push_back(term);
}

void 
ScoreFunction::operator+=(ScoreTermPtr term) {
	add_term(term);
}


VariedSpacerScoreFunction::VariedSpacerScoreFunction() {}

VariedSpacerScoreFunction::VariedSpacerScoreFunction(
		vector<string> spacers): my_spacers(spacers) {}

double
VariedSpacerScoreFunction::evaluate(
		ConstructConstPtr sgrna,
		EvaluatedScoreFunction &table) const {

	double score = 0;
	int N = my_spacers.size();

	// Return 0 if no spacers where specified.  Maybe this should throw, but my 
	// instinct is that it's better to return a sensible value.
	if(N == 0) { return 0; }
	
	// Score the sgRNA with all the spacers this score function knows about.
	table.clear();
	for(string spacer: my_spacers) {
		ConstructPtr sgrna_i = sgrna->copy();
		EvaluatedScoreFunction table_i;

		sgrna_i->domain("spacer")->seq(spacer);
		score += ScoreFunction::evaluate(sgrna_i, table_i);

		// Copy the score table entries into the "real" score table.
		for(auto &row: table_i) {
			row.name = spacer + ": " + row.name;
			row.term /= N;
			table.push_back(row);
		}
	}

	// Normalize the score by the number of spacers.
	return score / N;
}

vector<string>
VariedSpacerScoreFunction::spacers() const {
	return my_spacers;
}

void
VariedSpacerScoreFunction::spacers(vector<string> spacers) {
	my_spacers = spacers;
}

void
VariedSpacerScoreFunction::add_spacer(string spacer) {
	my_spacers.push_back(spacer);
}

void
VariedSpacerScoreFunction::operator+=(string spacer) {
	add_spacer(spacer);
}


ScoreTerm::ScoreTerm(string name, double weight):
	my_name(name), my_weight(weight) {}

string
ScoreTerm::name() const {
	return my_name;
}

void
ScoreTerm::name(string name) {
	my_name = name;
}

double
ScoreTerm::weight() const {
	return my_weight;
}

void
ScoreTerm::weight(double weight) {
	my_weight = weight;
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


FavorWildtypeTerm::FavorWildtypeTerm(
		ConstructConstPtr wt,
		vector<string> selection,
		double weight):

	ScoreTerm("favor_wt", weight),
	my_wt(wt),
	my_selection(selection) {}

double
FavorWildtypeTerm::evaluate(
		ConstructConstPtr sgrna,
		RnaFold const &apo_fold,
		RnaFold const &holo_fold) const {

	int wt_nucs = 0;
	int num_positions = 0;

	if(my_selection.empty()) {
		return 0;
	}

	for(string k: my_selection) {
		for(int i = 0; i < my_wt->domain(k)->len(); i++) {
			wt_nucs += (sgrna->domain(k)->seq()[i] == my_wt->domain(k)->seq()[i]);
			num_positions += 1;
		}
	}

	return static_cast<double>(wt_nucs) / num_positions;
}


ActiveMacrostateTerm::ActiveMacrostateTerm(double weight):
	ScoreTerm("active_macrostate", weight) {}

double
ActiveMacrostateTerm::evaluate(
		ConstructConstPtr sgrna,
		RnaFold const &apo_fold,
		RnaFold const &holo_fold) const {

	double p_apo = apo_fold.macrostate_prob(sgrna->active());
	double p_holo = holo_fold.macrostate_prob(sgrna->active());

	// Return the natural logarithm of the calculated probability so this score 
	// term can be properly added to others.
	return log(1 - p_apo) + log(p_holo);
}


// Most of these score terms have the same general form:
//
// Sum[
//  Prob{desired fold} × log(Enrichment{desired fold over undesired fold})
// ]
//
// A fold is a base pair for the "paired" score terms or an unpaired nucleotide 
// for the "unpaired" score terms.  For the "conditional" score terms, the 
// undesired fold is the same fold as the desired fold, but in the opposite 
// condition (e.g. apo vs. holo).  For the "always" score terms, the undesired 
// fold is anything that's not the desired fold, in the same condition.  So 
// Prob{undesired fold} = 1 - Prob{desired fold}.
//
// You can think of the probability term as counting how many nucleotides have 
// the desired fold.  You can think of the enrichment term as weighting that by 
// how much more likely the desired fold is compared to the undesired fold.  
// The logarithm is important because it allows the enrichments to be summed.
//
// In some cases, negative terms for undesired folds are included.  However, 
// these can be confusing to look at, because you have to remember that the 
// sign of the logarithm changes depending on which fold is more probable.  
// Written fully out, these negative terms would look like:
//
// + Prob{desired fold} × | log(Enrichment{desired over undesired}) |
// - Prob{undesired fold} × | log(Enrichment{undesired over desired}) |
//
// where |x| represents the absolute value of x.  It turns out that the above 
// expression has a more succinct form:
//
// (Prob{desired fold} + Prob{undesired fold})
//   × log(Enrichment{desired fold over undesired fold})

LigandSensitivityTerm::LigandSensitivityTerm(
		string name,
		vector<string> selection,
		double weight):

	ScoreTerm(name, weight),
	my_selection(selection) {}

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

	if(indices.empty()) {
		return 0;
	}

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


ConditionallyPairedTerm::ConditionallyPairedTerm(
		string name,
		ConditionEnum condition,
		vector<string> selection,
		vector<string> targets,
		double weight):

	ScoreTerm(name, weight),
	my_condition(condition),
	my_selection(selection),
	my_targets(targets) {}

double
ConditionallyPairedTerm::evaluate(
		ConstructConstPtr sgrna,
		RnaFold const &apo_fold,
		RnaFold const &holo_fold) const {

	vector<int> selection_indices = domains_to_indices(sgrna, my_selection);
	vector<int> target_indices = domains_to_indices(sgrna, my_targets);

	if(selection_indices.empty()) {
		return 0;
	}

	double sensitivity = 0;
	double p_apo, p_holo;

	for(int i: selection_indices) {
		for(int j: target_indices) {
			p_apo = apo_fold.base_pair_prob(i, j);
			p_holo = holo_fold.base_pair_prob(i, j);

			// The calculation below obscures the meaning of the score term.  I think 
			// the meaning is easier to see if you keep the following (mathematically 
			// equivalent) expression in mind:
			//
			// + Prob{paired, right condition} × log(Enrichment{right condition})
			// - Prob{paired, wrong condition} × log(Enrichment{wrong condition}),
			//
			// where Enrichment is how much more likely the nucleotide is to be 
			// paired with the relevant partner in the given condition. The logarithm 
			// is important because it allows us to sum these sub-expressions.

			if(p_apo > 0 and p_holo > 0) {
				switch(my_condition) {
					case ConditionEnum::APO:
						sensitivity += (p_apo + p_holo) * log(p_apo / p_holo);
						break;
					case ConditionEnum::HOLO:
						sensitivity += (p_holo + p_apo) * log(p_holo / p_apo);
						break;
				}
			}
		}
	}

	return sensitivity / selection_indices.size();
}


ConditionallyUnpairedTerm::ConditionallyUnpairedTerm(
		string name,
		ConditionEnum condition,
		vector<string> selection,
		double weight):

	ScoreTerm(name, weight),
	my_condition(condition),
	my_selection(selection) {}

double
ConditionallyUnpairedTerm::evaluate(
		ConstructConstPtr sgrna,
		RnaFold const &apo_fold,
		RnaFold const &holo_fold) const {

	vector<int> selection_indices = domains_to_indices(sgrna, my_selection);

	if(selection_indices.empty()) {
		return 0;
	}

	double sensitivity = 0;

	for(int i: selection_indices) {
		double p_apo = 1;   // Probability that nucleotide ``i`` is unpaired
		double p_holo = 1;  // in the given condition.

		for(int j = 0; j < sgrna->len(); j++) {
			p_apo -= apo_fold.base_pair_prob(i, j);
			p_holo -= holo_fold.base_pair_prob(i, j);
		}

		// The calculation below obscures the meaning of the score term.  I think 
		// the meaning is easier to see if you keep the following (mathematically 
		// equivalent) expression in mind:
		//
		// + Prob{unpaired, right condition} × log(Enrichment{right condition})
		// - Prob{unpaired, wrong condition} × log(Enrichment{wrong condition}),
		//
		// where Enrichment is how much more likely the nucleotide is to be 
		// unpaired in the given condition. The logarithm is important because it 
		// allows us to sum these sub-expressions.

		if(p_apo > 0 and p_holo > 0) {
			switch(my_condition) {
				case ConditionEnum::APO:
					sensitivity += (p_apo + p_holo) * log(p_apo / p_holo);
					break;
				case ConditionEnum::HOLO:
					sensitivity += (p_holo + p_apo) * log(p_holo / p_apo);
					break;
			}
		}
	}

	return sensitivity / selection_indices.size();
}


AlwaysPairedTerm::AlwaysPairedTerm(
		string name,
		vector<string> selection,
		vector<string> targets,
		double weight):

	ScoreTerm(name, weight),
	my_selection(selection),
	my_targets(targets) {}

double
AlwaysPairedTerm::evaluate(
		ConstructConstPtr sgrna,
		RnaFold const &apo_fold,
		RnaFold const &holo_fold) const {

	vector<int> selection_indices = domains_to_indices(sgrna, my_selection);
	vector<int> target_indices = domains_to_indices(sgrna, my_targets);

	if(selection_indices.empty()) {
		return 0;
	}

	double sensitivity = 0;

	for(int i: selection_indices) {
		double p_apo = 0;    // Probability that nucleotide ``i`` is paired
		double p_holo = 0;   // with any of its allowed partners.

		for(int j: target_indices) {
			p_apo += apo_fold.base_pair_prob(i, j);
			p_holo += holo_fold.base_pair_prob(i, j);
		}

		// Here, the probability of the undesired fold doesn't depend on the 
		// condition, it's just one minus the probability of the desired fold.  
		// This simplifies the basic score term expression even more, because the 
		// probabilities sum to one and the only remaining term is the enrichment 
		// ratio.  The enrichment ratios for both conditions are included.

		if(p_apo > 0 and p_apo < 1 and p_holo > 0 and p_holo < 1) {
			sensitivity += log(p_apo / (1 - p_apo));
			sensitivity += log(p_holo / (1 - p_holo));
		}
	}

	// Divide by two because there are two positive terms associated with each 
	// base pair, one for correctly forming in the apo condition and another for 
	// the holo condition.

	return sensitivity / selection_indices.size() / 2;
}


AlwaysUnpairedTerm::AlwaysUnpairedTerm(
		string name,
		vector<string> selection,
		double weight):

	ScoreTerm(name, weight),
	my_selection(selection) {}

double
AlwaysUnpairedTerm::evaluate(
		ConstructConstPtr sgrna,
		RnaFold const &apo_fold,
		RnaFold const &holo_fold) const {

	vector<int> selection_indices = domains_to_indices(sgrna, my_selection);

	if(selection_indices.empty()) {
		return 0;
	}

	double sensitivity = 0;

	for(int i: selection_indices) {
		double p_apo = 1;    // Probability that nucleotide ``i`` is unpaired
		double p_holo = 1;   // in the given condition.

		for(int j = 0; j < sgrna->len(); j++) {
			p_apo -= apo_fold.base_pair_prob(i, j);
			p_holo -= holo_fold.base_pair_prob(i, j);
		}

		// Here, the probability of the undesired fold doesn't depend on the 
		// condition, it's just one minus the probability of the desired fold.  
		// This simplifies the basic score term expression even more, because the 
		// probabilities sum to one and the only remaining term is the enrichment 
		// ratio.  The enrichment ratios for both conditions are included.

		if(p_apo > 0 and p_apo < 1 and p_holo > 0 and p_holo < 1) {
			sensitivity += log(p_apo / (1 - p_apo));
			sensitivity += log(p_holo / (1 - p_holo));
		}
	}

	// Divide by two because there are two positive terms associated with each 
	// nucleotide, one for being unpaired in the apo condition and another for 
	// the holo condition.

	return sensitivity / selection_indices.size() / 2;
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

