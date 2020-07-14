/*
Mass Prefix Tree.

Tree structure that mimics prefix trees for characters but uses floating
point numbers instead.

Original code from:
https://www.geeksforgeeks.org/trie-insert-and-search/

Modified by:
Zachary McGrath
*/


// Returns new trie node (initialized to NULLs)
struct TrieNode *getNode(void);

// If not present, inserts key into trie
// If the key is prefix of trie node, just marks leaf node
void insert(struct TrieNode *root, const char *key);

// Returns true if key presents in trie, else false
bool search(struct TrieNode *root, const char *key);
