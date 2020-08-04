#include "catch.hpp"
#include "../src/MassDawgNode.hpp"

TEST_CASE("Mass Dawg Node test cases"){
    MassDawgNode * mdn;
    MassDawgNode * child;

    float delta = 0.0001;

    SECTION("Initiating a Mass Dawg Node with no kmers gives a node with no kmers or children"){
        mdn = new MassDawgNode();

        REQUIRE(mdn->kmers.size() == 0);
        REQUIRE(mdn->children.size() == 0);
    }

    mdn = new MassDawgNode(100.1, 200.2, "ABC");
    SECTION("Creating a Mass Dawg Node with a kmer and 2 masses gives a node with a single kmer, 2 masses and no childre"){

        REQUIRE(mdn->kmers.size() == 1);
        REQUIRE(mdn->kmers[0] == "ABC");
        REQUIRE((double)abs(mdn->singlyMass - 100.1) < delta);
        REQUIRE((double)abs(mdn->doublyMass - 200.2) < delta);

        REQUIRE(mdn->children.size() == 0);
    }

    mdn->addKmer("DEF");
    SECTION("Adding a kmer to an existing node increases the size of the kmers vector and adds the kmer"){

        REQUIRE(mdn->kmers.size() == 2);
        REQUIRE(mdn->kmers[1] == "DEF");
        REQUIRE(mdn->kmers[0] == "ABC");
    }

    child = mdn->addChild(300.3, 400.4, "XYZ");
    SECTION("Adding a child to the node returns a child with the 2 masses and kmer passed in"){

        REQUIRE(mdn->children.size() == 1);
        REQUIRE(mdn->children[0] == child);

        REQUIRE(child->kmers.size() == 1);
        REQUIRE(child->children.size() == 0);
        REQUIRE((double)abs(child->singlyMass - 300.3) < delta);
        REQUIRE((double)abs(child->doublyMass - 400.4) < delta);
        REQUIRE(child->kmers[0] == "XYZ");
    }

    SECTION("Hashing of nodes return different values for nodes with different masses"){
        REQUIRE_NOTHROW(mdn->hash());
        REQUIRE_NOTHROW(child->hash());

        REQUIRE(mdn->hash() != child->hash());
    }

    delete mdn;
    delete child;
    mdn = new MassDawgNode(100.1, 200.2, "ABC");
    child = new MassDawgNode(300.3, 400.4, "XYZ");
    SECTION("Adding child by pointer does not throw an error and child can be found"){
        REQUIRE_NOTHROW(mdn->addChildByPointer(child));
        REQUIRE(mdn->children[0]->kmers[0] == "XYZ");
    }

    SECTION("Adding a child by pointer then trying to add it again will not add it the second time"){
        REQUIRE_NOTHROW(mdn->addChildByPointer(child));
        REQUIRE_NOTHROW(mdn->addChildByPointer(child));
        REQUIRE(mdn->children.size() == 1);
    }

    SECTION("Adding a child by pointer then adding another one with masses within the delta of .001 will not add the second child"){
        REQUIRE_NOTHROW(mdn->addChildByPointer(child));
        
        MassDawgNode * otherChild = new MassDawgNode(300.3001, 400.4001, "XYZ");
        REQUIRE_NOTHROW(mdn->addChildByPointer(otherChild));
        REQUIRE(mdn->children.size() == 1);
    }

    delete mdn;
    delete child;

}