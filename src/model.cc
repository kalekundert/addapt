#include <algorithm>

#include "model.hh"
#include "utils.hh"

namespace addapt {

Device::Device(string seq):
	my_seq(seq), my_context(make_shared<Context>()) {}

int
Device::len() const {
	return seq().length();
}

string
Device::seq() const {
	return my_context->before() + my_seq + my_context->after();
}

char
Device::seq(int index) const {
	index = normalize_index(seq(), index, IndexEnum::ITEM);
	return seq()[index];
}

int
Device::raw_len() const {
	return my_seq.length();
}

string
Device::raw_seq() const {
	return my_seq;
}

char
Device::raw_seq(int index) const {
	index = normalize_index(my_seq, index, IndexEnum::ITEM);
	return my_seq[index];
}

string
Device::macrostate(string name) const {
	if(my_macrostates.find(name) == my_macrostates.end()) {
		throw (f("no macrostate '%s'") % name).str();
	}
	string before(my_context->before().length(), '.');
	string after(my_context->after().length(), '.');
	return before + my_macrostates.at(name) + after;
}

Device::macrostate_view
Device::macrostates() const {
	return macrostate_view(*this);
}

void
Device::add_macrostate(string name, string macrostate) {
	if(macrostate.length() != my_seq.length()) {
		throw "constraint length doesn't match sequence length";
	}
	my_macrostates[name] = macrostate;
}

ContextConstPtr
Device::context() const {
	return my_context;
}

void
Device::context(ContextConstPtr context) {
	my_context = context;
}

void
Device::remove_context() {
	my_context = make_shared<Context>();
}

void
Device::mutate(int context_indep_index, char const mutation) {
	int index = normalize_index(my_seq, context_indep_index, IndexEnum::ITEM);
	my_seq[index] = mutation;
}

DevicePtr
Device::copy() const {
	DevicePtr other = std::make_shared<Device>(my_seq);
	other->my_macrostates = my_macrostates;
	other->my_context = my_context;
	return other;
}

void
Device::assign(DevicePtr other) {
	my_seq = other->my_seq;
	my_macrostates = other->my_macrostates;
	my_context = other->my_context;
}


Aptamer::Aptamer(string seq, string fold, double affinity):
	my_seq(seq), my_fold(fold), my_affinity(affinity) {}

string
Aptamer::seq() const {
	return my_seq;
}

string
Aptamer::fold() const {
	return my_fold;
}

double
Aptamer::affinity() const {
	return my_affinity;
}


Context::Context(string before, string after):
	my_before(before), my_after(after) {}

string
Context::before() const {
	return my_before;
}

string
Context::after() const {
	return my_after;
}


}
