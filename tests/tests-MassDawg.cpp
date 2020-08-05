#include <vector>

#include "catch.hpp"
#include "../src/MassDawg.hpp"

using namespace std;

bool hasString(vector<string> listOfString, string searching){
    for (string s: listOfString){
        if (s == searching) return true;
    }
    return false;
}


TEST_CASE("Testing Mass Dawg"){
    MassDawg * md;

    md = new MassDawg();

    string searchString1 = "ABCD";
    vector<float> singlySearchSeq1 = {200.2, 400.4, 600.6, 800.8};
    vector<float> doublySearchSeq1 = {100.1, 200.2, 300.3, 400.4};

    string searchString2 = "ABYZ";
    vector<float> singlySearchSeq2 = {200.2, 400.4, 700.7, 900.9};
    vector<float> doublySearchSeq2 = {100.1, 200.2, 350.35, 450.45};

    string searchString3 = "WXYZ";
    vector<float> singlySearchSeq3 = {100.1, 340.34, 700.7, 900.9};
    vector<float> doublySearchSeq3 = {50.05, 170.17, 350.35, 450.45};

    string searchString4 = "MNXY";
    vector<float> singlySearchSeq4 = {180.18, 260.26, 700.7, 900.9};
    vector<float> doublySearchSeq4 = {90.09, 130.13, 350.35, 450.45};

    string searchString5 = "DEFG";
    vector<float> singlySearchSeq5 = {150.15, 380.38, 750.75, 980.98};
    vector<float> doublySearchSeq5 = {75.075, 190.19, 325.325, 490.49};

    SECTION("One insert into graph throws no error and kmer is found in search"){
        REQUIRE_NOTHROW(md->insert(singlySearchSeq1, doublySearchSeq1, searchString1));
    }

    SECTION("Two insertions into graph throws no error and both kmers can be found"){
        REQUIRE_NOTHROW(md->insert(singlySearchSeq1, doublySearchSeq1, searchString1));
        REQUIRE_NOTHROW(md->insert(singlySearchSeq2, doublySearchSeq2, searchString2));

        REQUIRE(hasString(md->fuzzySearch(singlySearchSeq1, 0, 10), searchString1));
        REQUIRE(hasString(md->fuzzySearch(singlySearchSeq2, 0, 10), searchString2));

        REQUIRE(hasString(md->search(singlySearchSeq1, 10), searchString1));
        REQUIRE(hasString(md->search(singlySearchSeq2, 10), searchString2));
    }

    SECTION("Two insertions out of order does not throw exception"){
        REQUIRE_NOTHROW(md->insert(singlySearchSeq2, doublySearchSeq2, searchString2));
        REQUIRE_NOTHROW(md->insert(singlySearchSeq1, doublySearchSeq1, searchString1));
        md->finish();

        REQUIRE(hasString(md->search(singlySearchSeq1, 10), searchString1));
        REQUIRE(hasString(md->search(singlySearchSeq2, 10), searchString2));
    }

    SECTION("Two insertions with a common prefix and a search only until the common prefix returns both kmers"){
        REQUIRE_NOTHROW(md->insert(singlySearchSeq1, doublySearchSeq1, searchString1));
        REQUIRE_NOTHROW(md->insert(singlySearchSeq2, doublySearchSeq2, searchString2));

        vector<string> results = md->fuzzySearch({200.2, 400.4}, 0, 10);
        
        REQUIRE(hasString(results, "AB"));
        REQUIRE(hasString(md->search({200.2, 400.4}, 10), "AB"));
    }

    SECTION("Two insertions with no common prefix but a common suffix and subsequent search for it returns both kmers"){
        REQUIRE_NOTHROW(md->insert(singlySearchSeq3, doublySearchSeq3, searchString3));
        REQUIRE_NOTHROW(md->insert(singlySearchSeq2, doublySearchSeq2, searchString2));
        REQUIRE_NOTHROW(md->finish());

        vector<string> results = md->fuzzySearch(singlySearchSeq2, 0, 10);

        REQUIRE(hasString(results, searchString2));
        REQUIRE(hasString(results, searchString3));
    }

    SECTION("Two insertions with no common prefix but a common suffix and subsequent fuzzy search with a gap for it returns both kmers"){
        REQUIRE_NOTHROW(md->insert(singlySearchSeq3, doublySearchSeq3, searchString3));
        REQUIRE_NOTHROW(md->insert(singlySearchSeq2, doublySearchSeq2, searchString2));
        REQUIRE_NOTHROW(md->finish());

        vector<string> results = md->fuzzySearch({200.2, 700.7, 900.9}, 1, 10);

        REQUIRE(hasString(results, searchString2));
        REQUIRE(hasString(results, searchString3));
    }

    SECTION("Three insertions with two shared common suffixes and subsequent search returns those shared suffixes"){
        REQUIRE_NOTHROW(md->insert(singlySearchSeq5, doublySearchSeq5, searchString5));
        REQUIRE_NOTHROW(md->insert(singlySearchSeq1, doublySearchSeq1, searchString1));
        REQUIRE_NOTHROW(md->insert(singlySearchSeq2, doublySearchSeq2, searchString2));
        REQUIRE_NOTHROW(md->finish());

        vector<string> results = md->fuzzySearch({200.2, 400.4}, 0, 10);
        
        REQUIRE(hasString(results, "AB"));
    }

    SECTION("Three insertions with all three shared suffixes and subseqeunt search returns all three kmers"){
        REQUIRE_NOTHROW(md->insert(singlySearchSeq3, doublySearchSeq3, searchString3));
        REQUIRE_NOTHROW(md->insert(singlySearchSeq4, doublySearchSeq4, searchString4));
        REQUIRE_NOTHROW(md->insert(singlySearchSeq2, doublySearchSeq2, searchString2));
        REQUIRE_NOTHROW(md->finish());

        vector<string> results = md->fuzzySearch(singlySearchSeq3, 0, 10);
        
        REQUIRE(hasString(results, searchString3));
        REQUIRE(hasString(results, searchString4));
        REQUIRE(hasString(results, searchString2));
    }

    SECTION("Three insertions with no shared suffixes and full length search does not return all 3 full kmers"){
        REQUIRE_NOTHROW(md->insert(singlySearchSeq3, doublySearchSeq3, searchString3));
        REQUIRE_NOTHROW(md->insert(singlySearchSeq5, doublySearchSeq5, searchString5));
        REQUIRE_NOTHROW(md->insert(singlySearchSeq1, doublySearchSeq1, searchString1));
        REQUIRE_NOTHROW(md->finish());

        vector<string> results3 = md->fuzzySearch(singlySearchSeq3, 1, 10);
        vector<string> results5 = md->fuzzySearch(singlySearchSeq5, 1, 10);
        vector<string> results1 = md->fuzzySearch(singlySearchSeq1, 1, 10);
        
        // make sure that the other seqeuences are not in the other results
        REQUIRE(hasString(results3, searchString3));
        REQUIRE_FALSE(hasString(results3, searchString5));
        REQUIRE_FALSE(hasString(results3, searchString1));

        REQUIRE_FALSE(hasString(results5, searchString3));
        REQUIRE(hasString(results5, searchString5));
        REQUIRE_FALSE(hasString(results5, searchString1));

        REQUIRE_FALSE(hasString(results1, searchString3));
        REQUIRE_FALSE(hasString(results1, searchString5));
        REQUIRE(hasString(results1, searchString1));
    }
}