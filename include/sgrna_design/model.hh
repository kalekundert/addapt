#pragma once

#include <iostream>
#include <memory>
#include <vector>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/mem_fun.hpp>

#include <sgrna_design/utils.hh>

namespace sgrna_design {

class Domain;
using DomainPtr = std::shared_ptr<Domain>;
using DomainConstPtr = std::shared_ptr<Domain const>;
using DomainList = std::vector<DomainPtr>;

class Construct;
using ConstructPtr = std::shared_ptr<Construct>;
using ConstructConstPtr = std::shared_ptr<Construct const>;

struct Nucleotide {

	/// @brief The index of this nucleotide in the given construct.
	string domain;

	/// @brief The offset of this nucleotide relative to the start or end  
	/// (positive or negative, respectively) of the domain.
	int offset;

};

class Sequence {

public:

	/// @brief Return the sequence represented by this object.
	virtual string seq() const = 0;

	/// @brief Return the length of the sequence.
	int len() const;

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

	/// @brief Return the nucleotide at the given index of this domain.
	char seq(int) const;

	/// @brief Set the sequence of this domain.
	void seq(string const);

	/// @brief Make a point mutation in this domain.
	void mutate(int, char);

	/// @brief Insert a new sequence into this domain.
	/// @details The insertion is made immediately after the given index.
	void insert(int, char);

	/// @brief Insert a new sequence into this domain.
	/// @details The insertion is made immediately after the given index.
	void insert(int, string);

	/// @brief Delete the given index of this domain.
	void remove(int);

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

struct hash {};
using DomainMultiIndex = boost::multi_index_container<
	DomainPtr,
	boost::multi_index::indexed_by<
		boost::multi_index::random_access<>, // insertion order
		boost::multi_index::hashed_unique<   // lookup by name
			boost::multi_index::tag<hash>,
			boost::multi_index::const_mem_fun<Domain, string, &Domain::name> > > >;

class Construct : public Sequence {

public:

	/// @brief Default constructor.
	Construct();

	/// @brief Return a deep-copy of this construct.
	ConstructPtr copy() const;

	/// @brief Return the sequence of this construct.
	string seq() const;

	/// @brief Return the domains that make up this construct.
	DomainMultiIndex domains() const;

	/// @brief Return the domain with the given name.
	DomainPtr domain(string) const;

	/// @brief Return the domain with the given name.
	DomainPtr operator[](string) const;

	/// @brief Return the domain with the given name.
	DomainPtr operator[](DomainConstPtr) const;

	/// @brief Add the given domain to this construct.
	void add_domain(DomainPtr);

	/// @brief Add the given domain to this construct.
	void operator+=(DomainPtr);

	/// @brief Return the absolute index of the given nucleotide, which specifies 
	/// a domain and an index within that domain.
	int index(Nucleotide) const;

	/// @brief Return the index for the 5' side of the given domain.
	int index_5(string) const;
	int index_5(DomainConstPtr) const;

	/// @brief Return the index for the 3' side of the given domain.
	int index_3(string) const;
	int index_3(DomainConstPtr) const;

private:

	DomainMultiIndex my_domains;

};


std::ostream &operator<<(std::ostream &, Construct const &);
std::ostream &operator<<(std::ostream &, Domain const &);


}
