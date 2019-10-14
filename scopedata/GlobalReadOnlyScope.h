#include<unordered_map>

/**
 * Global scope of the VarDict. Contains configuration that must be available from all the classes and methods
 * and current run mode of VarDict for starting pipelines.
 * Must be initialized only once. Clear method created only for testing purposes.
 */
class GlobalReadOnlyScope {

    private:
		volatile static GlobalReadOnlyScope instance;
    	volatile static AbstractMode mode;

    public :
		static GlobalReadOnlyScope instance() {
        	return instance;
    	};

    	static synchronized void init(Configuration conf, unordered_map<string, int> chrLengths, string sample, string samplem,
                                         string ampliconBasedCalling, unordered_map<string, int> adaptorForward,
                                         unordered_map<string, int> adaptorReverse) {
        if (instance != null) {
            throw new IllegalStateException("GlobalReadOnlyScope was already initialized. Must be initialized only once.");
        }
        instance = new GlobalReadOnlyScope(conf, chrLengths, sample, samplem, ampliconBasedCalling, adaptorForward,
                adaptorReverse);
    	};


    	static AbstractMode getMode() {
        	return mode;
    	};

    	static synchronized void setMode(AbstractMode runMode) {
        	if (mode != null) {
            	throw new IllegalStateException("Mode was already initialized for GlobalReadOnlyScope. Must be initialized only once.");
        	}
        	mode = runMode;
    	};	

    /**
     * TEST usage only
     */
    	static synchronized void clear(){
        	instance = null;
        	mode = null;
    	};
		final Configuration conf;
  		final unordered_map<string, int> chrLengths;
  		final string sample;
  		final string samplem;
  		final string ampliconBasedCalling;
  		final PrinterType printerTypeOut;
  		final unordered_map<string, int> adaptorForward;
  		final unordered_map<string, int> adaptorReverse;

   		GlobalReadOnlyScope(Configuration conf, unordered_map<string, int> chrLengths, string sample, string samplem,
                               string ampliconBasedCalling, unordered_map<string, int> adaptorForward,
                               unordered_map<string, int> adaptorReverse) {
        	this->conf = conf;
        	this->chrLengths = chrLengths;
        	this->sample = sample;
        	this->samplem = samplem;
        	this->ampliconBasedCalling = ampliconBasedCalling;
        	this->printerTypeOut = conf.printerType;
        	this->adaptorForward = adaptorForward;
        	this->adaptorReverse = adaptorReverse;
    	};
}
