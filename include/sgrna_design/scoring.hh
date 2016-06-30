#pragma once

#include <iostream>
#include <memory>
#include <vector>

extern "C" {
  #include <ViennaRNA/data_structures.h>
}

#include <sgrna_design/model.hh>

namespace sgrna_design {

class ScoreFunction;
using ScoreFunctionPtr = std::shared_ptr<ScoreFunction>;

class ScoreTerm;
using ScoreTermPtr = std::shared_ptr<ScoreTerm>;
using ScoreTermList = std::vector<ScoreTermPtr>;

enum class LigandEnum {
	NONE,
	THEO,
};

enum class ConditionEnum {
	APO,
	HOLO,
};


class ScoreFunction {

public:

	/// @brief Default constructor.
	ScoreFunction();

	/// @brief Add a term to this score function.
	void operator+=(ScoreTermPtr);

	/// @brief Calculate a score for the given construct.
	double evaluate(ConstructPtr) const;

private:
	ScoreTermList my_terms;

};

/// @brief The interface to RNA secondary structure predictions.
class RnaFold {

public:

	/// @brief Return the probability that these two nucleotides will base pair 
	/// with each other.
	virtual double base_pair_prob(int, int) const = 0;

};

class ViennaRnaFold : public RnaFold {

public:

	/// @brief Predict how the construct will fold.
	ViennaRnaFold(ConstructConstPtr, LigandEnum=LigandEnum::NONE);

	/// @brief Free the ViennaRNA data structures.
	~ViennaRnaFold();

	/// @brief Return a string depicting the predicted structure for this 
	/// sequence.
	const char *base_pair_string() const;

	/// @brief Return the probability that these two nucleotides will base pair 
	/// with each other.
	double base_pair_prob(int, int) const;

private:

	ConstructConstPtr my_sgrna;
	// We need to store our own copy of the sgRNA sequence to prevent the pointer 
	// returned by c_str() from getting deleted.
	string my_seq;
	vrna_fold_compound_t *my_fc;
	char *my_fold;
};


class ScoreTerm {

public:

	/// @brief Default constructor.
	ScoreTerm();

	/// @brief Initialize the score term with a weight.
	ScoreTerm(double);

	/// @brief Return this score term's weight.
	double weight() const;

	/// @brief Set this score term's weight.
	void weight(double);

	/// @brief Calculate an unweighted value for this score term.
	virtual double evaluate(
			ConstructConstPtr, RnaFold const &, RnaFold const &) const = 0;

private:

	double my_weight;

};

class BasePairingTerm : public ScoreTerm {

public:

	BasePairingTerm(
			ConditionEnum, vector<string>, vector<string>, double=1.0);

	double evaluate(
			ConstructConstPtr, RnaFold const &, RnaFold const &) const;

private:

	ConditionEnum my_condition;
	vector<string> my_selection_a;
	vector<string> my_selection_b;

};

class LigandSensitivityTerm : public ScoreTerm {

public:

	LigandSensitivityTerm(
			vector<string>, double=1.0);

	double evaluate(
			ConstructConstPtr, RnaFold const &, RnaFold const &) const;

private:

	vector<string> my_selection;

};

class SpecificLigandSensitivityTerm : public ScoreTerm {

public:

	SpecificLigandSensitivityTerm(
			ConditionEnum, vector<string>, vector<string>, double=1.0);

	double evaluate(
			ConstructConstPtr, RnaFold const &, RnaFold const &) const;

private:

	ConditionEnum my_condition;
	vector<string> my_selection;
	vector<string> my_targets;


};

}

namespace std {

ostream&
operator<<(ostream&, const sgrna_design::ConditionEnum&);

}
