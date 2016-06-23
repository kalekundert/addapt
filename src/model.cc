#include <boost/format.hpp>
using boost::format;

#include "model.hh"
#include "utils.hh"

namespace sgrna {

char const *
Sequence::c_str() const {
	return seq().c_str();
}

int
Sequence::len() const {
	return seq().length();
}


Construct::Construct() {}

ConstructPtr
Construct::copy() const {
	// Make a new construct.
	ConstructPtr construct = std::make_shared<Construct>();

	// Fill it with copies of my domains.
	for(DomainPtr domain: my_domains) {
		construct += domain->copy();
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

DomainList
Construct::domains() const {
	return my_domains;
}

void
Construct::domain(DomainPtr domain) {
	my_domains.push_back(domain);
}

int
Construct::index(Nucleotide nuc) const {
	int index = 0;

	for(auto domain : my_domains) {
		if (nuc.domain == domain->name()) {
			return index + nuc.index;
		}
		index += domain->len();
	}

	throw format("No nucleotide named %s") % nuc.domain;
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

void
Domain::seq(string const seq) {
	my_seq = seq;
}

void
Domain::mutate(int index, char mutation) {
	my_seq[index] = mutation;
}

void
Domain::insert(int index, string insert) {
	my_seq.insert(index + 1, insert);
}

void
Domain::remove(int start, int end) {
	my_seq.erase(start, end - start);
}

void
Domain::replace(int start, int end, string insert) {
	my_seq.replace(start, end - start, insert);
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


void operator+=(ConstructPtr construct, DomainPtr domain) {
	construct->domain(domain);
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
