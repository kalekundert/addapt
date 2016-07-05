#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/mem_fun.hpp>

extern "C" {
  #include <ViennaRNA/data_structures.h>
}

#include <sgrna_design/model.hh>

namespace sgrna_design {

class ScoreFunction;
using ScoreFunctionPtr = std::shared_ptr<ScoreFunction>;

 struct EvaluatedScoreTerm { string name; double weight, term; };
using EvaluatedScoreFunction = std::vector<EvaluatedScoreTerm>;

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


class ScoreFunction {

public:

	/// @brief Default constructor.
	ScoreFunction();

	/// @brief Add a term to this score function.
	void add_term(ScoreTermPtr);

	/// @brief Add a term to this score function.
	void operator+=(ScoreTermPtr);

	/// @brief Calculate a score for the given construct.
	double evaluate(ConstructConstPtr) const;

	/// @brief Calculate a score for the given construct and fill in a table 
	/// containing the name, weight, and value of each score term.
	double evaluate(ConstructConstPtr, EvaluatedScoreFunction &) const;

private:
	ScoreTermList my_terms;

};

class ScoreTerm {

public:

	/// @brief Optionally initialize the score term with a name and a weight.
	ScoreTerm(string="", double=1.0);

	/// @brief Return this score term's name.
	string name() const;

	/// @brief Set this score term's name.
	void name(string);

	/// @brief Return this score term's weight.
	double weight() const;

	/// @brief Set this score term's weight.
	void weight(double);

	/// @brief Calculate an unweighted value for this score term.
	virtual double evaluate(
			ConstructConstPtr, RnaFold const &, RnaFold const &) const = 0;

private:

	string my_name;
	double my_weight;

};

class FavorWildtypeTerm : public ScoreTerm {

public:

	FavorWildtypeTerm(
			ConstructConstPtr, vector<string>, double=1.0);

	/// @brief The fraction of the nucleotides in the selected domains that have 
	/// been mutated.
	double evaluate(
			ConstructConstPtr, RnaFold const &, RnaFold const &) const;

private:

	ConstructConstPtr my_wt;
	vector<string> my_selection;

};

class LigandSensitivityTerm : public ScoreTerm {

public:

	LigandSensitivityTerm(
			string, vector<string>, double=1.0);

	double evaluate(
			ConstructConstPtr, RnaFold const &, RnaFold const &) const;

private:

	vector<string> my_selection;

};

class ConditionallyPairedTerm : public ScoreTerm {

public:

	ConditionallyPairedTerm(
			string, ConditionEnum, vector<string>, vector<string>, double=1.0);

	double evaluate(
			ConstructConstPtr, RnaFold const &, RnaFold const &) const;

private:

	ConditionEnum my_condition;
	vector<string> my_selection;
	vector<string> my_targets;

};

class ConditionallyUnpairedTerm : public ScoreTerm {

public:

	ConditionallyUnpairedTerm(
			string, ConditionEnum, vector<string>, double=1.0);

	double evaluate(
			ConstructConstPtr, RnaFold const &, RnaFold const &) const;

private:

	ConditionEnum my_condition;
	vector<string> my_selection;

};

class AlwaysPairedTerm : public ScoreTerm {

public:

	AlwaysPairedTerm(
			string, vector<string>, vector<string>, double=1.0);

	double evaluate(
			ConstructConstPtr, RnaFold const &, RnaFold const &) const;

private:

	ConditionEnum my_condition;
	vector<string> my_selection;
	vector<string> my_targets;

};

class AlwaysUnpairedTerm : public ScoreTerm {

public:

	AlwaysUnpairedTerm(
			string, vector<string>, double=1.0);

	double evaluate(
			ConstructConstPtr, RnaFold const &, RnaFold const &) const;

private:

	ConditionEnum my_condition;
	vector<string> my_selection;

};


}

namespace std {

ostream&
operator<<(ostream&, const sgrna_design::ConditionEnum&);

}
