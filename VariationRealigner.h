#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <vector>
#include <iostream>

#include "parseCigar.h"
#include "patterns.h"
#include "Variation.h"
#include "VariationUtils.h"
#include "Configuration.h"
#include "data/BaseInsertion.h"
#include <map>
#include <sstream>
#include <unordered_map>
#include <regex>
#include "stdio.h"

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/hts_defs.h"
#include <sys/time.h>
using namespace std;


/**
 * The step of pipeline which try to realign variations: softclips, indels, long insertions.
 */
class VariationRealigner implements Module<VariationData, RealignedVariationData>  {

	private:
		unordered_map<int, List<Sclip>> SOFTP2SV;
	    unordered_map<int, VariationMap<string, Variation>> nonInsertionVariants;
	    unordered_map<int, VariationMap<string, Variation>> insertionVariants;
	    unordered_map<int, unordered_map<string, int>> positionToInsertionCount;
	    unordered_map<int, unordered_map<string, int>> positionToDeletionsCount;
	    unordered_map<int, int> refCoverage;
	    unordered_map<int, Sclip> softClips5End;
	    unordered_map<int, Sclip> softClips3End;
	    ReferenceResource referenceResource;
	    Reference reference;
	    Region region;
	    Set<String> splice;
	    string chr;
	    int maxReadLength;
	    double duprate;
	    string[] bams;
	    string bam;
	    unordered_map<int, unordered_map<string, int>> mnp;
	    SVStructures svStructures;
	    VariantPrinter variantPrinter;

    public:
		Scope<RealignedVariationData> process(Scope<VariationData> scope);

    	void initFromScope(Scope<VariationData> scope);
    	void filterAllSVStructures(); 
    	void filterSV(List<Sclip> svList_sva);
    	Cluster checkCluster(List<Mate> mates, int rlen) ;
    	void adjustMNP() ;
    	void realignIndels();
    	void realigndel(String[] bamsParameter, Map<Integer, Map<String, Integer>> positionToDeletionsCount); 
    	String realignins(Map<Integer, Map<String, Integer>> positionToInsertionCount) ;
    	void realignlgdel(List<Sclip> svfdel, List<Sclip> svrdel) ;
    	void realignlgins30() 
    	void realignlgins(List<Sclip> svfdup, List<Sclip> svrdup) ;
    	List<SortPositionDescription> fillAndSortTmp(Map<Integer, Map<String, Integer>> changes); 
};


class SortPositionDescription {
	public:
        int position;
        String descriptionString;
        int count;

        SortPositionDescription(int position, String descriptionString, int count); 
    }

public Match35 find35match(string seq5, string seq3);
boolean noPassingReads(string chr, int start, int end, string[] bams);
public boolean ismatch(string seq1, string seq2, int dir);
public boolean ismatch(string seq1, string seq2, int dir, int MM);
public static boolean islowcomplexseq(string seq);
static int count(string str, char chr);
public static BaseInsertion adjInsPos(int bi, string ins, unordered_map<int, char> ref);
BaseInsertion findbi(string seq, int position, unordered_map<int, char> ref, final int dir, string chr);
int findbp(string sequence, int startPosition, unordered_map<Integer, char> ref, int direction, string chr) ;
void adjRefCnt(Variation tv, Variation ref, int len);
void adjRefFactor(Variation ref, double factor_f);
void addVarFactor(Variation vref, double factor_f);
MismatchResult findMM5(unordered_map<int, char> ref, int position, string wupseq);
MismatchResult findMM3(unordered_map<int, char> ref, int p, string sanpseq);

class MismatchResult {
        private :
			final List<Mismatch> mismatches;
       		final List<Integer> scp;
        	final int nm;
        	final int misp;
        	final string misnt;

        public :
			MismatchResult(List<Mismatch> mismatches, List<Integer> scp, int nm, int misp, string misnt);
       		List<int> getScp() ;
        	List<Mismatch> getMismatches() ;
        	int getNm() ;
       		int getMisp() ;
        	string getMisnt() ;
};

class Mismatch {
		public:
       		string mismatchSequence; //$mm
       		int mismatchPosition; //$mp
       		int end; //# or 5
        	Mismatch(string mismatchSequence, int mismatchPosition, int end); 
       
};
public static boolean ismatchref(String sequence, unordered_map<int, char> ref, int position, int dir) ;
public static boolean ismatchref(string sequence, unordered_map<int, char> ref, int position, int dir, int MM); 
void rmCnt(Variation vref, Variation tv); 
