package org.broadinstitute.hellbender.tools.walkers.variantutils;

        import org.broadinstitute.hellbender.CommandLineProgramTest;
        import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
        import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
        import org.testng.annotations.Test;

        import java.util.Arrays;
        import java.util.Collections;

        import static org.testng.Assert.*;

/**
 * Created by gauthier on 10/2/17.
 */
public class ReblockGVCFIntegrationTest extends CommandLineProgramTest {
    //TODO: for the love of God do NOT commit the reference!  This is just to test that we can match the GATK3 results on chromosome 1
    String b37KGReference = getToolTestDataDir() + "human_g1k_v37.fasta";

    @Test
    public void testJustOneSample() throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-L 1:69485-69791 -O %s -R " + b37KGReference +
                        " -V " + getToolTestDataDir() + "gvcfForReblocking.g.vcf -RGQthreshold 20" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(getToolTestDataDir() + "testJustOneSample.expected.g.vcf"));
        spec.executeTest("testJustOneSample", this);
    }

    @Test
    //covers non-ref GT correction, but not non-ref AD when non-ref is not called
    public void testProductionGVCF() throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-L 1:1-1000000 -O %s -R " + b37KGReference +
                        " -V " + getToolTestDataDir() + "NA12878.prod.chr1snippet.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false" +
                        " -RGQ-threshold 40",
                Arrays.asList(getToolTestDataDir() + "testProductionGVCF.expected.g.vcf"));
        spec.executeTest("testProductionGVCF", this);
    }

    @Test
    public void testOneSampleDropLows() throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-dropLowQuals -L 1:69485-69791 -O %s -R " + b37KGReference +
                        " -V " + getToolTestDataDir() + "gvcfForReblocking.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(getToolTestDataDir() + "testOneSampleDropLows.expected.g.vcf"));
        spec.executeTest("testOneSampleDropLows", this);
    }

    @Test
    public void testNonRefADCorrection() throws Exception {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                "-O %s -R " + b38_reference_20_21 +
                        " -V " + getToolTestDataDir() + "nonRefAD.g.vcf" +
                        " --" + StandardArgumentDefinitions.ADD_OUTPUT_VCF_COMMANDLINE + " false",
                Arrays.asList(getToolTestDataDir() + "testNonRefADCorrection.expected.g.vcf"));
        spec.executeTest("testNonRefADCorrection", this);
    }
}