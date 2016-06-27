#include <algorithm>

#include <boost/format.hpp>
using boost::format;

#include <sgrna_design/model.hh>
#include <sgrna_design/utils.hh>

namespace sgrna_design {

int
Sequence::len() const {
	return seq().length();
}

int
Sequence::normalize_index(int index, bool between) const {
	int normalized_index = index;

	// If the user gave a negative index, interpret it as counting backward from 
	// the end of the sequence.
	if(index < 0) {
		normalized_index += len() + between;
	}

	// Make sure the index refers to a position that actually exists in the 
	// sequence.  The maximum index is one greater if the index refers to the 
	// positions between the nucleotides rather than the nucleotides themselves.  
	if(normalized_index < 0 or normalized_index > len() - (between? 0:1)) {
		throw (f("no index '%d' in '%s'") % index % seq()).str();
	}

	// Return the normalized index.
	return normalized_index;
}

pair<int,int>
Sequence::normalize_range(int start, int end) const {
	// Resolve negative and out-of-bounds indices.
	start = normalize_index(start, true);
	end = normalize_index(end, true);

	// Work out which index is lower and which is higher.
	int normalized_start = std::min(start, end);
	int normalized_end = std::max(start, end);

	// Return the normalized indices.
	return {normalized_start, normalized_end};
}


Domain::Domain(
		string const name,
		string const seq,
		ColorEnum color,
		StyleEnum style):
	
	my_name(name),
	my_seq(seq),
	my_color(color),
	my_style(style) {}

DomainPtr
Domain::copy() const {
	return std::make_shared<Domain>(my_name, my_seq, my_color, my_style);
}

string
Domain::name() const {
	return my_name;
}

void
Domain::name(string const name) {
	my_name = name;
}

string
Domain::seq() const {
	return my_seq;
}

char
Domain::seq(int index) const {
	return my_seq[index];
}

void
Domain::seq(string const seq) {
	my_seq = seq;
}

void
Domain::mutate(int index, char mutation) {
	index = normalize_index(index);
	my_seq[index] = mutation;
}

void
Domain::insert(int index, string insert) {
	index = normalize_index(index, true);
	my_seq.insert(index, insert);
}

void
Domain::remove(int start, int end) {
	auto indices = normalize_range(start, end);
	my_seq.erase(indices.first, indices.second - indices.first);
}

void
Domain::replace(int start, int end, string insert) {
	auto indices = normalize_range(start, end);
	my_seq.replace(indices.first, indices.second - indices.first, insert);
}

ColorEnum
Domain::color() const {
	return my_color;
}

void
Domain::color(ColorEnum color) {
	my_color = color;
}

StyleEnum
Domain::style() const {
	return my_style;
}

void
Domain::style(StyleEnum style) {
	my_style = style;
}


Construct::Construct() {}

ConstructPtr
Construct::copy() const {
	// Make a new construct.
	ConstructPtr construct = std::make_shared<Construct>();

	// Fill it with copies of my domains.
	for(DomainPtr domain: my_domains) {
		*construct += domain->copy();
	}

	// Return the new construct.
	return construct;
}

string
Construct::seq() const {
	string seq;

	for(DomainPtr domain: my_domains) {
		seq += domain->seq();
	}

	return seq;
}

DomainMultiIndex
Construct::domains() const {
	return my_domains;
}

DomainPtr
Construct::domain(string name) const {
	auto it = my_domains.get<hash>().find(name);
	if (it == my_domains.get<hash>().end()) {
		throw (f("no such domain '%s'") % name).str();
	}
	return *it;
}

DomainPtr
Construct::operator[](string name) const {
	return domain(name);
}

DomainPtr
Construct::operator[](DomainConstPtr domain) const {
	return this->domain(domain->name());
}

void
Construct::add_domain(DomainPtr domain) {
	my_domains.push_back(domain);
}

void
Construct::operator+=(DomainPtr domain) {
	add_domain(domain);
}

int
Construct::index(Nucleotide nuc) const {
	int index = 0;

	for(DomainConstPtr domain : my_domains) {
		if (nuc.domain == domain->name()) {
			return index + nuc.offset;
		}
		index += domain->len();
	}

	throw (format("No nucleotide named %s") % nuc.domain).str();
}

int
Construct::index_5(string name) const {
	int index = 0;

	for(DomainConstPtr domain : my_domains) {
		if (name == domain->name()) {
			return index;
		}
		index += domain->len();
	}

	throw (format("No domain named %s") % name).str();
}

int
Construct::index_5(DomainConstPtr domain) const {
	return index_5(domain->name());
}

int
Construct::index_3(string name) const {
	return index_3(domain(name));
}

int
Construct::index_3(DomainConstPtr domain) const {
	return index_5(domain) + domain->len() - 1;
}


std::ostream &operator<<(std::ostream &out, Construct const &construct) {
	for(DomainPtr domain : construct.domains()) {
		out << *domain;
	}
	return out;
}

std::ostream &operator<<(std::ostream &out, Domain const &domain) {
	out << color(domain.seq(), domain.color(), domain.style());
	return out;
}


}
