#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include "utils.hh"

namespace addapt {

class Construct;
using ConstructPtr = std::shared_ptr<Construct>;
using ConstructConstPtr = std::shared_ptr<Construct const>;

class Aptamer;
using AptamerPtr = std::shared_ptr<Aptamer>;
using AptamerConstPtr = std::shared_ptr<Aptamer const>;

class Construct {

public:

	/// @brief Default constructor.
	Construct(string);

	/// @brief Return the sequence of this construct.
	string seq() const;

	/// @brief Return the nucleotide at the given position of this construct.
	char seq(int) const;

	/// @brief Return constraints that define the given macrostate.
	string macrostate(string) const;

	/// @brief Return all the macrostates associated with this construct.
	unordered_map<string,string> macrostates() const;

	/// @brief Add constraints that define a particular macrostate.
	void add_macrostate(string, string);

	/// @brief Return the length of this construct.
	int len() const;

	/// @brief Make a point mutation in this construct.
	void mutate(int, char const);

	/// @brief Return a deep-copy of this construct.
	ConstructPtr copy() const;

	/// @brief Make this construct equivalent to the given one.
	void assign(ConstructPtr);

private:

	string my_seq;
	unordered_map<string,string> my_macrostates;

};

class Aptamer {

public:

	/// @brief Construct with a sequence, fold, and ΔG (kcal/mol).
	Aptamer(string, string, double);

	/// @brief Return the sequence of this aptamer.
	string seq() const;

	/// @brief Return a pseudo-dot-bracket string specifying the base pairs 
	/// formed by this aptamer in the holo condition.
	string fold() const;

	/// @brief Return the affinity of this aptamer for its ligand (i.e. its  
	/// dissociation constant in units of μM).
	double affinity() const;

private:

	string my_seq;
	string my_fold;
	double my_affinity;

};


}
