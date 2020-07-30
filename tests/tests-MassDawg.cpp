#include "catch.hpp"
#include "../src/MassDawg.hpp"

TEST_CASE("Testing Mass Dawg"){
    MassDawg * md;

    md = new MassDawg();

    SECTION("One insert into graph throws no error and kmer is found in search"){

    }

    SECTION("Two insertions in order into graph throws no error and both kmers can be found"){

    }

    SECTION("Two insertions out of order throws exception"){

    }

    SECTION("Two insertions with a common prefix and a search only until the common prefix returns both kmers"){

    }

    SECTION("Two insertions with no common prefix but a common suffix and subsequent search for it returns both kmers"){

    }

    SECTION("Three insertions in order with two shared common suffixes and subsequent search returns those shared suffixes"){

    }

    SECTION("Three insertions in order with all three shared suffixes and subseqeunt search returns all three kmers"){

    }

    SECTION("Three insertions in order with no shared suffixes and full length search does not return all 3 full kmers"){
        
    }
}