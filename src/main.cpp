#include <stdio.h>
#include <vector>
#include <iostream>

#include "MassDawg.hpp"

using namespace std;

int main(){
    // create a mass dawg
    MassDawg * md = new MassDawg();

    // add a sequence
    vector<float> doubly {100.1, 200.2, 300.3};
    vector<float> singly {200.2, 400.4, 600.6};
    string sequence = "ABC";

    // insert this sequcne
    cout << "inserting ABC into graph";
    md->insert(singly, doubly, sequence);

    // show it before the next addition
    md->show();

    doubly = {100.1, 200.2, 300.3, 400.4};
    singly = {200.2, 400.4, 600.6, 800.8};
    sequence = "ABCD";

    cout << "inserting ABCD into graph";
    md->insert(singly, doubly, sequence);

    // show the graph
    md->show();

    doubly = {100.1, 200.4, 300.6, 400.4};
    singly = {200.2, 400.4, 600.6, 800.8};
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

    cout << "\nsearching for a sequence with gaps\n";
    searching = {300.3, 800.8};
    results = md->fuzzySearch(searching, 2, 10);
    cout << "Number of results: " + to_string(results.size()) + "\n";

    for (string i: results){
        cout << i + "\n";
    }
    searching = {100.1, 400.4};
    results = md->search(searching, 10);
    for (string i: results) cout << i + "\n";

    delete md;
}