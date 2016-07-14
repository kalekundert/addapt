#pragma once

#include <iostream>
#include <memory>
#include <vector>
#include <list>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/mem_fun.hpp>

extern "C" {
  #include <ViennaRNA/data_structures.h>
}

#include "model.hh"

namespace addapt {

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

enum class FavorableEnum {
	NO,
	YES,
};


/// @brief The interface to RNA secondary structure predictions.
class RnaFold {

public:

	/// @brief Return the probability that these two nucleotides will base pair 
	/// with each other.
	virtual double base_pair_prob(int, int) const = 0;

	/// @brief Return the probability that the device passed to the 
	/// deviceor will fold into the given macrostate, defined by a hard 
	/// constraint string.
	virtual double macrostate_prob(string) const = 0;

};

class ViennaRnaFold : public RnaFold {

public:

	/// @brief Predict how the device will fold.
	ViennaRnaFold(DeviceConstPtr, AptamerConstPtr=nullptr);

	/// @brief Free the ViennaRNA data structures.
	~ViennaRnaFold();

	/// @brief Return the probability that these two nucleotides will base pair 
	/// with each other.
	double base_pair_prob(int, int) const;

	/// @brief Return the probability that the device passed to the 
	/// deviceor will fold into the given macrostate, defined by a hard 
	/// constraint string.
	double macrostate_prob(string) const;

private:

	vrna_fold_compound_t *make_fold_compound(bool) const;

private:

	DeviceConstPtr my_sgrna;
	AptamerConstPtr my_aptamer;

	// We (used to) need to store our own copy of the sgRNA sequence to ensure 
	// that the pointer returned by c_str() is valid.
	string my_seq;

	// We need to store a list of all the fold compound data structures we end up 
	// creating so we can deallocate them all.
	mutable list<vrna_fold_compound_t *> my_fcs;

	// We need a fold compound object to cache the base-pair probability matrix.
	mutable vrna_fold_compound_t *my_bppm_fc;
};


class ScoreFunction {

public:

	/// @brief Default deviceor.
	ScoreFunction();

	/// @brief Calculate a score for the given device.
	double evaluate(DeviceConstPtr) const;

	/// @brief Calculate a score for the given device and fill in a table 
	/// containing the name, weight, and value of each score term.
	virtual double evaluate(DeviceConstPtr, EvaluatedScoreFunction &) const;

	/// @brief Add a term to this score function.
	void add_term(ScoreTermPtr);

	/// @brief Add a term to this score function.
	void operator+=(ScoreTermPtr);

	/// @brief Return the aptamer being used by this score function.
	AptamerConstPtr aptamer() const;

	/// @brief Set the aptamer being used by this score function.
	void aptamer(AptamerConstPtr);

private:
	ScoreTermList my_terms;
	AptamerConstPtr my_aptamer;

};

/// @brief Calculate a score in the context of several different spacer 
/// sequences and return the average.
///
/// @details The purpose of this score function is to find sequences that are 
/// somehow robust to changes in the spacer sequence.
class VariedSpacerScoreFunction : public ScoreFunction {

public:

	/// @brief Default deviceor.
	VariedSpacerScoreFunction();

	/// @brief Deviceor that accepts a list of spacer sequences.
	VariedSpacerScoreFunction(vector<string>);

	/// @brief Calculate a score for the given device and fill in a table 
	/// containing the name, weight, and value of each score term.
	virtual double evaluate(DeviceConstPtr, EvaluatedScoreFunction &) const;

	/// @brief Return the spacer sequences being used by this score function.
	vector<string> spacers() const;

	/// @brief Set the spacer sequences that will be used by this score function.
	void spacers(vector<string>);

	/// @brief Add a spacer sequence to those being considered by this score 
	/// function.
	void add_spacer(string);

	/// @brief Add a spacer sequence to those being considered by this score 
	/// function.
	void operator+=(string);

private:

	/// @brief A list of spacer sequences.  The sequences in this list should 
	/// each be 20 nucleotides long, but this is not checked.
	vector<string> my_spacers;

};


class ScoreTerm {

public:

	/// @brief Optionally initialize the score term with a name and a weight.
	ScoreTerm(string="", double=1.0);

	/// @brief Calculate an unweighted value for this score term.
	virtual double evaluate(
			DeviceConstPtr, RnaFold const &, RnaFold const &) const = 0;

	/// @brief Return this score term's name.
	string name() const;

	/// @brief Set this score term's name.
	void name(string);

	/// @brief Return this score term's weight.
	double weight() const;

	/// @brief Set this score term's weight.
	void weight(double);

private:

	string my_name;
	double my_weight;

};

class MacrostateProbTerm : public ScoreTerm {

public:

	/// @brief Initialize the score term with a fold, a condition, and an 
	/// indication of whether or not we want that fold in that condition.
	MacrostateProbTerm(string, ConditionEnum, FavorableEnum=FavorableEnum::YES);

	/// @brief Calculate the log probability that the given device adopts the 
	/// given fold in the given condition.
	double evaluate(DeviceConstPtr, RnaFold const &, RnaFold const &) const;

private:
		string my_macrostate;
		ConditionEnum my_condition;
		FavorableEnum my_favorable;

};


}

namespace std {

ostream&
operator<<(ostream&, const addapt::ConditionEnum&);

ostream&
operator<<(ostream&, const addapt::FavorableEnum&);

}
