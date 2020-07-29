#include <stdio.h>
#include <vector>
#include <iostream>

#include "MassDawg.hpp"

using namespace std;

int main(){
    // create a mass dawg
    MassDawg * md = new MassDawg();

    // add a sequence
    vector<float> singly {100.1, 200.2, 300.3};
    vector<float> doubly {200.2, 400.4, 600.6};
    string sequence = "ABC";

    // insert this sequcne
    cout << "inserting ABC into graph";
    md->insert(singly, doubly, sequence);

    // show it before the next addition
    md->show();

    singly = {100.1, 200.2, 300.3, 400.4};
    doubly = {200.2, 400.4, 600.6, 800.8};
    sequence = "ABCD";

    cout << "inserting ABCD into graph";
    md->insert(singly, doubly, sequence);

    // show the graph
    md->show();

    singly = {100.1, 200.4, 300.6, 400.4};
    doubly = {200.2, 400.4, 600.6, 800.8};
    sequence = "AXYD";

    cout << "inserting AXYD int graph";
    md->insert(singly, doubly, sequence);

    md->finish();

    md->show();

    cout << "searching for a sequence with no gaps\n";
    vector<float> searching = {100.1, 400.4, 600.6, 800.8};
    vector<string> results = md->fuzzySearch(searching, 0, 10);
    cout << "Number of results: " + to_string(results.size()) + "\n";

    for (int i = 0; i < results.size(); i++){
        cout << results[i] + "\n"; 
    }

}