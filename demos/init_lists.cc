#include <iostream>
#include <vector>
#include <memory>
#include <initializer_list>

using namespace std;

void f1(vector<string> vec) {
    for(auto x: vec) {
        cout << x << ' ';
    }
		cout << endl;
}

struct C1 {
    C1(vector<string> vec) {
        for(auto x : vec) {
            cout << x << ' ';
        }
				cout << endl;
    }
};

struct C2 {
    C2(initializer_list<string> init) {
        for(auto x : init) {
            cout << x << ' ';
        }
				cout << endl;
    }
};


int main(int argc, const char * argv[]) {
    f1({"foo", "bar", "baz"});
    C1({"foo", "bar", "baz"});
		C2({"foo", "bar", "baz"});
		auto p = make_shared<C2>({"foo", "bar", "baz"});
    //auto ptr = make_shared<Func>({"foo", "bar", "baz"}); // won't compile.
    //auto ptr = make_shared<Func>{"foo", "bar", "baz"}; // nor this.
    return 0;
}
