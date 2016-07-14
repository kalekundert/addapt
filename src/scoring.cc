#include <algorithm>
#include <cmath>

#include <boost/algorithm/string.hpp>

extern "C" {
  #include <ViennaRNA/structure_utils.h>
  #include <ViennaRNA/part_func.h>
  #include <ViennaRNA/fold.h>
	#include <ViennaRNA/constraints.h>
}

#include "scoring.hh"

namespace addapt {

ViennaRnaFold::ViennaRnaFold(DeviceConstPtr device, AptamerConstPtr aptamer):

	my_device(device),
	my_aptamer(aptamer),
	my_seq(device->seq()),
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
ViennaRnaFold::macrostate_prob(string constraint) const {
	vrna_fold_compound_t *fc = make_fold_compound(false);

	// Calculate the free energy for the whole ensemble.
	double g_tot = vrna_pf(fc, NULL);

	// Add a constraint that defines the given macrostate.
	vrna_constraints_add(fc, constraint.c_str(),
			VRNA_CONSTRAINT_DB_DEFAULT | VRNA_CONSTRAINT_DB_ENFORCE_BP);

	// Calculate the free energy for the given macrostate.
	double g_active = vrna_pf(fc, NULL);

	// Return the probability that the device will be in the given macrostate 
	// at equilibrium.
	double kT = fc->exp_params->kT / 1000;
	return exp((g_tot - g_active) / kT);
}

vrna_fold_compound_t *
ViennaRnaFold::make_fold_compound(bool compute_bppm) const {
	// Make sure the device hasn't changed since this engine was created.
	assert(my_device->len() == my_seq.length());

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

	// Add the aptamer, if we were given one.
	if (my_aptamer) {
		double kT = fc->exp_params->kT / 1000;
		vrna_sc_add_hi_motif(
				fc,
				my_aptamer->seq().c_str(),
				my_aptamer->fold().c_str(),
				kT * log(my_aptamer->affinity() / 1e6),
				VRNA_OPTION_DEFAULT);
	}

	return fc;
}


ScoreFunction::ScoreFunction() {}

double
ScoreFunction::evaluate(DeviceConstPtr device) const {
	EvaluatedScoreFunction table;
	return evaluate(device, table);
}
		
double
ScoreFunction::evaluate(
		DeviceConstPtr device,
		EvaluatedScoreFunction &table) const {

	double score = 0;

	table.clear();

	if(my_contexts.empty()) {
		score += evaluate_terms(device, table);
	}
	else {
		DevicePtr scratch_device = device->copy();
		for(auto context: my_contexts) {
			scratch_device->context(context.second);
			score += evaluate_terms(scratch_device, table, context.first + ": ");
		}
	}

	return score;
}

double
ScoreFunction::evaluate_terms(
		DeviceConstPtr device,
		EvaluatedScoreFunction &table,
		string term_prefix) const {

	double score = 0;

	ViennaRnaFold apo_fold(device);
	ViennaRnaFold holo_fold(device, my_aptamer);

	for(ScoreTermPtr term: my_terms) {
		EvaluatedScoreTerm eval;
		eval.name = term_prefix + term->name();
		eval.weight = term->weight();
		eval.term = term->evaluate(device, apo_fold, holo_fold);
		table.push_back(eval);
		score += eval.weight * eval.term;
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

AptamerConstPtr
ScoreFunction::aptamer() const {
	return my_aptamer;
}

void
ScoreFunction::aptamer(AptamerConstPtr aptamer) {
	my_aptamer = aptamer;
}

ContextConstPtr
ScoreFunction::context(string name) const {
	return my_contexts.at(name);
}

void 
ScoreFunction::add_context(string name, ContextConstPtr context) {
	my_contexts[name] = context;
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


MacrostateProbTerm::MacrostateProbTerm(
		string macrostate,
		ConditionEnum condition,
		FavorableEnum favorable):

	ScoreTerm("active_macrostate"),
	my_macrostate(macrostate),
	my_condition(condition),
	my_favorable(favorable) {

	string desc;
	desc += (my_condition == ConditionEnum::APO)? "apo" : "holo";
	desc += ": ";
	desc += (my_favorable == FavorableEnum::NO)? "not " : "";
	desc += macrostate;
	name(desc);
}

double
MacrostateProbTerm::evaluate(
		DeviceConstPtr device,
		RnaFold const &apo_fold,
		RnaFold const &holo_fold) const {

	// Get a pointer to the right folding engine.
	RnaFold const *apropos_fold;
	switch(my_condition) {
		case ConditionEnum::APO: apropos_fold = &apo_fold; break;
		case ConditionEnum::HOLO: apropos_fold = &holo_fold; break;
	}

	// Calculate the probability of adopting this fold in this condition.
	string constraint = device->macrostate(my_macrostate);
	double macrostate_prob = apropos_fold->macrostate_prob(constraint);

	// Invert the probability if we want to avoid this fold in this condition.
	switch(my_favorable) {
		case FavorableEnum::YES: break;
		case FavorableEnum::NO: macrostate_prob = 1 - macrostate_prob; break;
	}

	// Return the natural logarithm of the calculated probability to improve 
	// dynamic range and to allow these terms to be summed.
	return log(macrostate_prob);
}


}

namespace std {

ostream &
operator<<(ostream& out, const addapt::ConditionEnum& condition) {
	switch(condition) {
		case addapt::ConditionEnum::APO: out << "APO"; break;
		case addapt::ConditionEnum::HOLO: out << "HOLO"; break;
	}
	return out;
}

ostream &
operator<<(ostream& out, const addapt::FavorableEnum& condition) {
	switch(condition) {
		case addapt::FavorableEnum::YES: out << "FAVORABLE"; break;
		case addapt::FavorableEnum::NO: out << "UNFAVORABLE"; break;
	}
	return out;
}

}

