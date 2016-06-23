#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include "utils.hh"

namespace sgrna {

class Construct;
using ConstructPtr = std::shared_ptr<Construct>;
using ConstructConstPtr = std::shared_ptr<Construct const>;

class Domain;
using DomainPtr = std::shared_ptr<Domain>;
using DomainConstPtr = std::shared_ptr<Domain const>;
using DomainList = std::vector<DomainPtr>;

struct Nucleotide {

	/// @brief The name of the domain this nucleotide is part of.
	string domain;

	/// @brief The index of this nucleotide in the given domain.
	/// @details This index may be either positive or negative.  If negative, it 
	/// is assumed to count from the end of the domain.
	int index;

};

using NucleotideList = std::vector<Nucleotide>;

class Sequence {

public:

	/// @brief Return the sequence represented by this object.
	virtual string seq() const = 0;

	/// @brief Return a C-style string (i.e. a pointer to an array of const 
	/// chars) containing the sequence represented by this object.  This is the 
	/// data type expected by all the ViennaRNA functions.
	char const * c_str() const;

	/// @brief Return the length of the sequence.
	int len() const;

};

class Construct : public Sequence {

public:

	/// @brief Default constructor.
	Construct();

	/// @brief Return a deep-copy of this construct.
	ConstructPtr copy() const;

	/// @brief Return the sequence of this construct.
	string seq() const;

	/// @brief Return the domains that make up this construct.
	DomainList domains() const;

	/// @brief Add the given domain to this construct.
	void domain(DomainPtr);

	/// @brief Return the absolute index of the given nucleotide, which specifies 
	/// a domain and an index within that domain.
	int index(Nucleotide) const;

private:

	//struct name {};
	//multi_index_container<
	//	DomainPtr,
	//	indexed_by<
	//		random_access<>, // insertion order
	//	  hashed_unique<   // lookup by name
	//			tag<name>,
	//			const_mem_fun<Domain, string, &Domain::name>
	//> my_domains;

	DomainList my_domains;

};

class Domain : public Sequence {

public:

	/// @brief Default constructor.
	Domain(
			string const name,
			string const seq="",
			ColorEnum color=ColorEnum::NORMAL,
			StyleEnum style=StyleEnum::NORMAL);

	/// @brief Return a deep-copy of this domain.
	DomainPtr copy() const;

	/// @brief Return the name of this domain.
	string name() const;

	/// @brief Set the name of this domain.
	void name(string const);

	/// @brief Return the sequence of this domain.
	string seq() const;

	/// @brief Set the sequence of this domain.
	void seq(string const);

	/// @brief Make a point mutation in this domain.
	void mutate(int, char);

	/// @brief Insert a new sequence into this domain.
	/// @details The insertion is made immediately after the given index.
	void insert(int, string);

	/// @brief Delete a section of this domain.
	void remove(int, int);

	/// @brief Replace a section of this domain with a new sequence.
	void replace(int, int, string);

public:

	/// @brief Return the color of this domain.
	ColorEnum color() const;

	/// @brief Set the color of this domain.
	void color(ColorEnum);

	/// @brief Return the test style of this domain.
	StyleEnum style() const;

	/// @brief Set the test style of this domain.
	void style(StyleEnum);

private:

	string my_name;
	string my_seq;
	ColorEnum my_color;
	StyleEnum my_style;

};

void operator+=(ConstructPtr, DomainPtr);
std::ostream &operator<<(std::ostream &, Construct const &);
std::ostream &operator<<(std::ostream &, Domain const &);

}
