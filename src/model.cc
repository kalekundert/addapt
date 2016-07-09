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


Domain::Domain(
		string const name,
		ColorEnum color,
		StyleEnum style):
	
	Domain(name, "", "", color, style) {}

Domain::Domain(
		string const name,
		string const seq,
		ColorEnum color,
		StyleEnum style):
	
	Domain(name, seq, "", color, style) {}

Domain::Domain(
		string const name,
		string const seq,
		string const active,
		ColorEnum color,
		StyleEnum style):
	
	my_name(name),
	my_color(color),
	my_style(style) {
	
	Domain::seq(seq, active);
}

DomainPtr
Domain::copy() const {
	return std::make_shared<Domain>(
			my_name, my_seq, my_active, my_color, my_style);
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
Domain::seq(string const seq, string const active) {
	my_seq = seq;

	if(active.empty()) {
		my_active = std::string(seq.length(), '.');
	}
	else {
		my_active = active;
	}

	if(my_seq.length() != my_active.length()) {
		throw (f("'%s' and '%s' are different lengths.") % my_seq % my_active).str();
	}
}

string
Domain::active() const {
	return my_active;
}

void
Domain::mutate(int index, char mutation) {
	index = normalize_index(seq(), index, IndexEnum::ITEM);
	my_seq[index] = mutation;
}

void
Domain::insert(int index, char insertion) {
	insert(index, string(1, insertion));
}

void
Domain::insert(int index, string insertion) {
	index = normalize_index(seq(), index, IndexEnum::BETWEEN);
	my_seq.insert(index, insertion);
}

void
Domain::remove(int index) {
	index = normalize_index(seq(), index, IndexEnum::ITEM);
	my_seq.erase(index, 1);
};

void
Domain::remove(int start, int end) {
	auto indices = normalize_range(seq(), start, end, IndexEnum::BETWEEN);
	my_seq.erase(indices.first, indices.second - indices.first);
}

void
Domain::replace(int start, int end, string insert) {
	auto indices = normalize_range(seq(), start, end, IndexEnum::BETWEEN);
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

string
Construct::active() const {
	string cst;

	for(DomainPtr domain: my_domains) {
		cst += domain->active();
	}

	return cst;
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
