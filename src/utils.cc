#include "utils.hh"

namespace sgrna_design {

int
normalize_index(string sequence, int index, IndexEnum meaning) {
	int normalized_index = index;
	int const seq_len = sequence.length();

	// If the user gave a negative index, interpret it as counting backward from 
	// the end of the sequence.
	if(index < 0) {
		normalized_index += seq_len + (meaning == IndexEnum::BETWEEN);
	}

	// Make sure the index refers to a position that actually exists in the 
	// sequence.  The maximum index is one greater if the index refers to the 
	// positions between the nucleotides rather than the nucleotides themselves.  
	int const max_index = seq_len - (meaning == IndexEnum::ITEM);
	if(normalized_index < 0 or normalized_index > max_index) {
		throw (f("no index '%d' in '%s'") % index % sequence).str();
	}

	// Return the normalized index.
	return normalized_index;
}

pair<int,int>
normalize_range(string sequence, int start, int end, IndexEnum between) {
	// Resolve negative and out-of-bounds indices.
	start = normalize_index(sequence, start, between);
	end = normalize_index(sequence, end, between);

	// Work out which index is lower and which is higher.
	int normalized_start = std::min(start, end);
	int normalized_end = std::max(start, end);

	// Return the normalized indices.
	return {normalized_start, normalized_end};
}

string
color(
		string str,
		ColorEnum color=ColorEnum::NORMAL,
		StyleEnum style=StyleEnum::NORMAL) {

	auto fmt = boost::format("\033[%d;%dm%s\033[0;0m")
		% static_cast<int>(style)
		% static_cast<int>(color)
		% str;

	return fmt.str();
}

}
