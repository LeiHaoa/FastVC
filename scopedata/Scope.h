
/**
 * Common scope of data must be storing between steps of VarDict pipeline.
 * @param <T> data of current step of pipeline
 */
public class Scope<T> {

    public:
		final string bam;
    	final Region region;
    	final Reference regionRef;
    	final ReferenceResource referenceResource;
    	final int maxReadLength;
    	final set<string> splice;
    	final VariantPrinter out;

    	final T data;

    	Scope(string bam, Region region, Reference regionRef, ReferenceResource referenceResource, int maxReadLength,
    	             set<string> splice, VariantPrinter out, T data) {
    	    this.bam = bam;
    	    this.region = region;
    	    this.regionRef = regionRef;
    	    this.referenceResource = referenceResource;
    	    this.maxReadLength = maxReadLength;
    	    this.splice = splice;
    	    this.data = data;
    	    this.out = out;
    	}

    	public Scope(Scope<?> inheritableScope, T data) {
    	    this(inheritableScope.bam, inheritableScope.region, inheritableScope.regionRef, inheritableScope.referenceResource,
    	            inheritableScope.maxReadLength, inheritableScope.splice, inheritableScope.out, data);
    	}
}
