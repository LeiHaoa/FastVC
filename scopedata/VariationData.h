#include<string.h>
#include<unordered_map>
#include<set>
/**
 * The data created after CigarParser step in pipeline. Used for process Variation in
 * Variation Realigner and Structural Variants analysis (realign it and searching for structural variants).
 */
public class VariationData {
    public: 
		unordered_map<int, Variationunordered_map<string, Variation>> nonInsertionVariants;
    	unordered_map<int, Variationunordered_map<string, Variation>> insertionVariants;
    	unordered_map<int, unordered_map<string, int>> positionToInsertionCount;
    	unordered_map<int, unordered_map<string, int>> positionToDeletionsCount;
    	SVStructures svStructures;
    	unordered_map<int, int> refCoverage;
    	unordered_map<int, Sclip> softClips5End;
    	unordered_map<int, Sclip> softClips3End;
    	int maxReadLength;
    	set<string> splice;
    	unordered_map<int, unordered_map<string, int>> mnp;
    	unordered_map<string, int[]> spliceCount;
    	double duprate;

}

