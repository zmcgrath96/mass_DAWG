#include <stdio.h>
#include <vector>
#include <iostream>

#include "MassDawg.hpp"

using namespace std;

int main(){

    MassDawg * md = new MassDawg();

    vector<float> singly {
        132.04776143499998,
        203.084875435,
        306.09406043499996,
        363.11552443499994,
        476.19958843499995,
        575.268002435,
        646.305116435,
        733.3371444349999
    };

    vector<float> doubly {
        66.52751893499999,
        102.04607593499999,
        153.55066843499998,
        182.06140043499997,
        238.60343243499997,
        288.137639435,
        323.65619643499997,
        367.17221043499995
    };

    string searching = "MACGLVAS";

    md->insert(singly, doubly, searching);

    vector<float> allMasses {132.04776143499998,
        203.084875435,
        306.09406043499996,
        363.11552443499994,
        476.19958843499995,
        575.268002435,
        646.305116435,
        733.3371444349999,
        66.52751893499999,
        102.04607593499999,
        153.55066843499998,
        182.06140043499997,
        238.60343243499997,
        288.137639435,
        323.65619643499997,
        367.17221043499995
    };

    // get the results and print them
    vector<string> results = md->fuzzySearch(allMasses, 0, 20);

    for (string result: results){
        cout << result + "\n";
    }

    return 1;
}