#include <cmath>

#include "scoring.hh"

namespace sgrna {

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

	my_fc = vrna_fold_compound(my_seq.c_str(), NULL, VRNA_OPTION_PF);
	vrna_exp_params_rescale(my_fc, NULL);

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
		double aptamer_dg = my_fc->exp_params->kT * log(aptamer_kd_uM / 1e6);
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
RnaFold::base_pair_prob(Nucleotide a, Nucleotide b) const {
	int i = my_sgrna->index(a);
	int j = my_sgrna->index(b);
	return my_fc->exp_matrices->probs[my_fc->iindx[i] - j];
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


PairedBase::PairedBase(double weight, Nucleotide nuc)
	: ScoreTerm(weight), my_nuc(nuc) {}

Nucleotide
PairedBase::nuc() const {
	return my_nuc;
}

void
PairedBase::nuc(Nucleotide nuc) {
	my_nuc = nuc;
}

NucleotideList
PairedBase::no_lig_partners() const {
	return my_no_lig_partners;
}

void
PairedBase::no_lig_partner(Nucleotide nuc) {
	my_no_lig_partners.push_back(nuc);
}

NucleotideList
PairedBase::lig_partners() const {
	return my_lig_partners;
}

void
PairedBase::lig_partner(Nucleotide nuc) {
	my_lig_partners.push_back(nuc);
}

double
PairedBase::evaluate(
		ConstructConstPtr sgrna,
		RnaFold const &no_lig,
		RnaFold const &lig) const {
	return 0;
}

}

