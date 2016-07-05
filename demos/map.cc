// You can assign initializer_lists to map elements.

#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <vector>

using namespace std;

template <typename T>
ostream& operator<< (ostream& out, const vector<T>& vec) {
  if ( !vec.empty() ) {
    out << '[';
    copy(vec.begin(), vec.end(), ostream_iterator<T>(out, ", "));
    out << "\b\b]";
  }
  return out;
}


int main ()
{
	using ints = vector<int>;
  map<char,ints> my_map;

  my_map['a'] = {1,2,3};
	my_map['b'] = {4,5};
	my_map['c'] = {6};

  cout << "mymap['a'] is " << my_map['a'] << '\n';
  cout << "mymap['b'] is " << my_map['b'] << '\n';
  cout << "mymap['c'] is " << my_map['c'] << '\n';
  cout << "mymap['d'] is " << my_map['d'] << '\n';

  return 0;
}
