package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.*;
import org.fusesource.leveldbjni.All;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

import static org.testng.Assert.*;

public class ReblockGVCFUnitTest {
    private final ReblockGVCF reblocker = new ReblockGVCF();
    private final static Allele LONG_REF = Allele.create("ACTG", true);
    private final static Allele DELETION = Allele.create("A", false);

    @Test
    public void testCleanUpHighQualityVariant() {

    }

    @Test
    public void testReblockVariant() {

    }

    @Test
    public void testMakeGQ0RefCall() {
        final Genotype g = makeG("sample1", LONG_REF, DELETION, 11, 0, 37);
        final VariantContext toBeNoCalled = makeDeletionVC("lowQualVar", Arrays.asList(LONG_REF, DELETION), LONG_REF.length(), g);
        final VariantContext noCalled = reblocker.makeGQ0RefCall(toBeNoCalled)
    }

    //TODO: these are duplicated from PosteriorProbabilitiesUtilsUnitTest but PR #4947 modifies VariantContextTestUtils, so I'll do some refactoring before the second of the two is merged
    private Genotype makeG(final String sample, final Allele a1, final Allele a2, final int... pls) {
        return new GenotypeBuilder(sample, Arrays.asList(a1, a2)).PL(pls).make();
    }

    private VariantContext makeDeletionVC(final String source, final List<Allele> alleles, final int refLength, final Genotype... genotypes) {
        final int start = 10;
        final int stop = start+refLength-1;
        return new VariantContextBuilder(source, "1", start, stop, alleles).genotypes(Arrays.asList(genotypes)).filters((String)null).make();
    }

}