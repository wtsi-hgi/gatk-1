package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.*;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.DefaultGATKVariantAnnotationArgumentCollection;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKAnnotationArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc.FixedAFCalculatorProvider;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.HaplotypeCallerGenotypingEngine;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceMode;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.IndexedSampleList;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;
import org.broadinstitute.hellbender.utils.logging.OneShotLogger;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.writers.GVCFWriter;
import picard.cmdline.programgroups.OtherProgramGroup;
import picard.cmdline.programgroups.VariantManipulationProgramGroup;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import java.util.*;

/**
 * Perform joint genotyping on gVCF files produced by HaplotypeCaller
 *
 * <p>
 * GenotypeGVCFs merges gVCF records that were produced as part of the Best Practices workflow for variant discovery
 * (see Best Practices documentation for more details) using the '-ERC GVCF' or '-ERC BP_RESOLUTION' mode of the
 * HaplotypeCaller, or result from combining such gVCF files using CombineGVCFs. This tool performs the multi-sample
 * joint aggregation step and merges the records together in a sophisticated manner: at each position of the input
 * gVCFs, this tool will combine all spanning records, produce correct genotype likelihoods, re-genotype the newly
 * merged record, and then re-annotate it.</p>
 *
 * <h3>Input</h3>
 * <p>
 * One or more HaplotypeCaller gVCFs to genotype.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A combined, genotyped VCF.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk ReblockGVCF \
 *   -R reference.fasta \
 *   --variant sample1.g.vcf \
 *
 *   -O sample1.reblocked.g.vcf
 * </pre>
 *
 * <h3>Caveats</h3>
 * <p>Only single-sample gVCF files produced by HaplotypeCaller can be used as input for this tool.</p>
 * <p>By default this tool only passes through annotations used by VQSR.  A different set of annotations can be specified with the usual -A argument.  Annotation groups are not supported by this tool.</p>
 *
 * <h3>Special note on ploidy</h3>
 * <p>This tool assumes diploid genotypes.</p>
 *
 */
@ExperimentalFeature
@CommandLineProgramProperties(summary = "Compress a single-sample GVCF from HaplotypeCaller by merging homRef blocks using new GQ band parameters",
        oneLineSummary = "Condenses homRef blocks in a single-sample GVCF",
        programGroup = OtherProgramGroup.class,
        omitFromCommandLine = true)
@DocumentedFeature
public class ReblockGVCF extends VariantWalker {

    private static String GVCF_BLOCK = "GVCFBlock";
    private static int PLOIDY_TWO = 2;  //FIXME: hack to avoid extra arg collections

    private VariantContextWriter vcfWriter;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which variants should be written")
    private File outputFile;

    @Argument(fullName="includeNonVariantSites", shortName="allSites", doc="Include loci found to be non-variant after genotyping")
    public boolean INCLUDE_NON_VARIANTS = true;

    @ArgumentCollection
    public GenotypeCalculationArgumentCollection genotypeArgs = new GenotypeCalculationArgumentCollection();

    /**
     * Which annotations to copy over for the reblocked output gVCF file.
     */
    @Advanced
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to recompute.  The single value 'none' removes the default annotations")
    protected List<String> annotationsToUse = Arrays.asList(VCFConstants.DEPTH_KEY, VCFConstants.RMS_MAPPING_QUALITY_KEY, GATKVCFConstants.READ_POS_RANK_SUM_KEY, GATKVCFConstants.MAP_QUAL_RANK_SUM_KEY, GATKVCFConstants.RAW_RMS_MAPPING_QUALITY_KEY);

    @Advanced
    @Argument(fullName="GVCFGQBands", shortName="GQB", doc="Exclusive upper bounds for reference confidence GQ bands " +
            "(must be in [1, 100] and specified in increasing order)")
    public List<Integer> GVCFGQBands = new ArrayList<>();
    {
        GVCFGQBands.add(20); GVCFGQBands.add(100);
    };

    @Advanced
    @Argument(fullName="drop-low-qual-variants", shortName="drop-low-quals", doc="Exclude variants and homRef blocks that are GQ0 from the reblocked GVCF to save space")
    protected boolean dropLowQuals = false;

    @Advanced
    @Argument(fullName="RGQ-threshold-to-drop", shortName="RGQ-threshold", doc="Reference genotype quality (PL[0]) value below which variant sites will be converted to GQ0 homRef calls")
    protected double RGQthreshold = 0.0;

    @Advanced
    @Argument(fullName="do-qual-score-approximation", shortName="do-qual-approx", doc="Add necessary INFO field annotation to perform QUAL approximation downstream")
    protected boolean doQualApprox = false;

    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate. Note that dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    private HaplotypeCallerArgumentCollection hcArgs = new HaplotypeCallerArgumentCollection();
    // the genotyping engine
    private HaplotypeCallerGenotypingEngine genotypingEngine;
    // the INFO field annotation key names to remove
    private final List<String> infoFieldAnnotationKeyNamesToRemove = new ArrayList<>();

    // INFO Header names that require alt alleles
    final Set<String> infoHeaderAltAllelesLineNames = new LinkedHashSet<>();

    public void onTraversalStart() {
        VCFHeader inputHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> inputHeaders = inputHeader.getMetaDataInSortedOrder();

        final Set<VCFHeaderLine> headerLines = new HashSet<>(inputHeaders);
        headerLines.addAll(getDefaultToolVCFHeaderLines());

        // create the genotyping engine
        hcArgs.outputMode = OutputMode.EMIT_ALL_SITES;
        hcArgs.annotateAllSitesWithPLs = true;
        hcArgs.genotypeArgs = new GenotypeCalculationArgumentCollection(genotypeArgs);
        hcArgs.emitReferenceConfidence = ReferenceConfidenceMode.GVCF;   //this is important to force emission of all alleles at a multiallelic site
        genotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, new IndexedSampleList(inputHeader.getGenotypeSamples()), FixedAFCalculatorProvider.createThreadSafeProvider(hcArgs), true);

        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_QUAL_APPROX_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.VARIANT_DEPTH_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MAPPING_QUALITY_DEPTH));

        // Remove GCVFBlocks
        headerLines.removeIf(vcfHeaderLine -> vcfHeaderLine.getKey().startsWith(GVCF_BLOCK));

        headerLines.remove(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.DOWNSAMPLED_KEY));
        headerLines.remove(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.HAPLOTYPE_SCORE_KEY));
        headerLines.remove(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY));
        headerLines.remove(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerLines.remove(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));


        if ( INCLUDE_NON_VARIANTS ) {
            // Save INFO header names that require alt alleles
            for ( final VCFHeaderLine headerLine : headerLines ) {
                if (headerLine instanceof VCFInfoHeaderLine ) {
                    if (((VCFInfoHeaderLine) headerLine).getCountType() == VCFHeaderLineCount.A) {
                        infoHeaderAltAllelesLineNames.add(((VCFInfoHeaderLine) headerLine).getID());
                    }
                }
            }
        }

        if ( dbsnp.dbsnp != null  ) {
            VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);
        }

        VariantContextWriter writer = createVCFWriter(outputFile);

        try {
            vcfWriter = new GVCFWriter(writer, GVCFGQBands, genotypeArgs.samplePloidy);
        } catch ( IllegalArgumentException e ) {
            throw new CommandLineException.BadArgumentValue("GQBands", "are malformed: " + e.getMessage());
        }        vcfWriter.writeHeader(new VCFHeader(headerLines, inputHeader.getGenotypeSamples()));

        logger.info("Notice that the -ploidy parameter is ignored in " + getClass().getSimpleName() + " tool as this is automatically determined by the input variant files");
    }

    // get VariantContexts from input gVCFs and regenotype
    public void apply(VariantContext variant, ReadsContext reads, ReferenceContext ref, FeatureContext features) {
        final VariantContext newVC = regenotypeVC(ref, variant);
        if (newVC != null) {
            vcfWriter.add(newVC);
        }
    }

    /**
     * Re-genotype (and re-annotate) a VariantContext
     *
     * @param ref            the ref context
     * @param originalVC     the combined genomic VC
     * @return a new VariantContext or null if the site turned monomorphic and we don't want such sites
     */
    protected VariantContext regenotypeVC(final ReferenceContext ref, final VariantContext originalVC) {
        if ( originalVC == null ) {
            throw new IllegalArgumentException("originalVC cannot be null");
        } else if (!isProperlyPolymorphic(originalVC) && !INCLUDE_NON_VARIANTS) {
            return null;
        }

        VariantContext result = originalVC;

        //Pass back ref-conf homRef sites/blocks to be combined by the GVCFWriter
        if (isHomRefBlock(result)) {
            return processHomRefBlock(result);
        }

        //don't need to calculate quals for sites with no data whatsoever or sites already genotyped homRef
        if (result.getAttributeAsInt(VCFConstants.DEPTH_KEY,0) > 0 && !isHomRefCall(result)) {
            final GenotypeLikelihoodsCalculationModel model = result.getType() == VariantContext.Type.INDEL
                    ? GenotypeLikelihoodsCalculationModel.INDEL
                    : GenotypeLikelihoodsCalculationModel.SNP;
            result = genotypingEngine.calculateGenotypes(originalVC, model, null);
        }

        if (result == null || (!isProperlyPolymorphic(result) && !INCLUDE_NON_VARIANTS)) {
            return null;
        }

        //variants with PL[0] less than threshold get turned to homRef with PL=[0,0,0], shouldn't get INFO attributes
        //make sure we can call het variants with GQ >= RGQthreshold in joint calling downstream
        if(shouldBeReblocked(result)) {
            return reblockVariant(result, originalVC);
        }
        //high quality variant
        else {
            Map<String, Object> attrMap = new HashMap<>();
            Map<String, Object> origMap = originalVC.getAttributes();
            //pare down annotations to just VQSR features
            for(final String annotation : annotationsToUse) {
                if(origMap.containsKey(annotation))
                    attrMap.put(annotation, origMap.get(annotation));
            }
            if (doQualApprox && result.getGenotype(0).hasPL()) {
                attrMap.put(GATKVCFConstants.RAW_QUAL_APPROX_KEY, result.getGenotype(0).getPL()[0]);
                int varDP = QualByDepth.getDepth(result.getGenotypes(), null);
                if (varDP == 0) {  //prevent QD=Infinity case
                    varDP = result.getAttributeAsInt(VCFConstants.DEPTH_KEY, Integer.MAX_VALUE);
                }
                attrMap.put(GATKVCFConstants.VARIANT_DEPTH_KEY, varDP);
            }
            VariantContextBuilder builder = new VariantContextBuilder(result);
            builder.attributes(attrMap);

            boolean allelesNeedSubsetting = false;
            List<Allele> allelesToDrop = new ArrayList<>();
            if (dropLowQuals) {
                //only put in alleles that are called if we're dropping low quality variants (mostly because this can introduce GVCF gaps if deletion alleles are dropped)
                for (final Allele currAlt : result.getAlternateAlleles()) {
                    boolean foundMatch = false;
                    for (final Allele gtAllele : result.getGenotype(0).getAlleles()) {
                        if (gtAllele.equals(currAlt, false)) {
                            foundMatch = true;
                            break;
                        }
                        if (gtAllele.equals(Allele.NON_REF_ALLELE)) {
                            if (dropLowQuals) { //don't regenotype, just drop it -- this is a GQ 0 case if ever I saw one
                                return null;
                            } else {
                                GenotypeBuilder gb = makeGQ0RefCall(result, attrMap);
                                return builder.alleles(Arrays.asList(result.getReference(), Allele.NON_REF_ALLELE)).unfiltered().log10PError(VariantContext.NO_LOG10_PERROR).attributes(attrMap).genotypes(gb.make()).make();
                            }
                        }
                    }
                    if (!foundMatch && !currAlt.isSymbolic()) {
                        allelesNeedSubsetting = true;
                        allelesToDrop.add(currAlt);
                    }
                }
            }
            //remove any AD reads for the non-ref
            int nonRefInd = result.getAlleleIndex(Allele.NON_REF_ALLELE);
            boolean genotypesWereModified = false;
            final ArrayList<Genotype> genotypesArray = new ArrayList<>();
            GenotypesContext newGenotypes = result.getGenotypes();
            Genotype g = result.getGenotype(0);
            if(g.hasAD()) {
                int[] ad = g.getAD();
                if (ad.length >= nonRefInd && ad[nonRefInd] > 0) { //only initialize a builder if we have to
                    genotypesWereModified = true;
                    GenotypeBuilder gb = new GenotypeBuilder(g);
                    ad[nonRefInd] = 0;
                    gb.AD(ad).DP((int) MathUtils.sum(ad));
                    genotypesArray.add(gb.make());
                    newGenotypes = GenotypesContext.create(genotypesArray);
                }
            }
            else {
                genotypesArray.add(g);
            }

            //we're going to approximate MQ_DP with the site-level DP (should be informative and uninformative reads), which is pretty safe because it will only differ if reads are missing MQ
            attrMap.put(GATKVCFConstants.MAPPING_QUALITY_DEPTH, originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY,0));

            if(allelesNeedSubsetting) {
                List<Allele> newAlleleSet = new ArrayList<>();
                for(final Allele a : result.getAlleles()) {
                    newAlleleSet.add(Allele.create(a,false));
                }
                newAlleleSet.removeAll(allelesToDrop);
                builder.alleles(newAlleleSet);
                if(!genotypesWereModified) {
                    builder.genotypes(AlleleSubsettingUtils.subsetAlleles(result.getGenotypes(), PLOIDY_TWO, result.getAlleles(), newAlleleSet, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, result.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0)));
                }
                else {  //again, only initialize a builder if we have to
                    builder.genotypes(AlleleSubsettingUtils.subsetAlleles(newGenotypes, PLOIDY_TWO, result.getAlleles(), newAlleleSet, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN, result.getAttributeAsInt(VCFConstants.DEPTH_KEY, 0)));
                }
                //only trim if we're subsetting alleles, and we only subset if we're allowed to drop sites, as per the -drop-low-quals arg
                return GATKVariantContextUtils.reverseTrimAlleles(builder.attributes(attrMap).unfiltered().make());
            }
            return builder.attributes(attrMap).genotypes(newGenotypes).unfiltered().make();
        }
    }

    /**
     * Determines whether the provided VariantContext has real alternate alleles.
     *
     * There is a bit of a hack to handle the <NON-REF> case because it is not defined in htsjdk.Allele
     * We check for this as a biallelic symbolic allele.
     *
     * @param vc  the VariantContext to evaluate
     * @return true if it has proper alternate alleles, false otherwise
     */
    private boolean isProperlyPolymorphic(final VariantContext vc) {
        //obvious cases
        if (vc == null || vc.getAlternateAlleles().isEmpty()) {
            return false;
        } else if (vc.isBiallelic()) {
            return !(vc.getAlternateAllele(0).equals(Allele.SPAN_DEL) ||
                    vc.getAlternateAllele(0).equals(GATKVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED) ||
                    vc.isSymbolic());
        } else {
            return true;
        }
    }

    private boolean isHomRefBlock(final VariantContext result) {
        return result.getLog10PError() == VariantContext.NO_LOG10_PERROR;
    }

    //treat homRef "calls", i.e. annotated variants with non-symbolic alt alleles and homRef genotypes, differently from het/homVar calls or homRef blocks
    private boolean isHomRefCall(final VariantContext result) {
        return result.getGenotype(0).isHomRef() && result.getLog10PError() != VariantContext.NO_LOG10_PERROR;
    }

    private VariantContext processHomRefBlock(final VariantContext result) {
        if (dropLowQuals && (result.getGenotype(0).getGQ() <= RGQthreshold || result.getGenotype(0).getGQ() == 0)) {
            return null;
        }
        else if (result.getGenotype(0).isCalled() && result.getGenotype(0).isHomRef()) {
            return result;
        }
        else if (!result.getGenotype(0).isCalled() && result.getGenotype(0).hasPL() && result.getGenotype(0).getPL()[0] == 0) {
            return result;
        }
        else
            return null;
    }

    private boolean shouldBeReblocked(final VariantContext result) {
        return result.getGenotype(0).getPL()[0] < RGQthreshold || result.getGenotype(0).isHomRef();
    }

    //"reblock" a variant by converting its genotyping to homRef, changing PLs, adding reblock END tags and other attributes
    private VariantContext reblockVariant(final VariantContext result, final VariantContext originalVC) {
        if(dropLowQuals && !isHomRefCall(result)) {
            return null;
        }
        Map<String, Object> attrMap = new HashMap<>();
        Allele newRef = result.getReference();
        Genotype g = result.getGenotype(0);
        GenotypeBuilder gb = new GenotypeBuilder(g);
        //NB: If we're dropping a deletion allele, then we need to trim the reference and add an END tag with the vc stop position
        if(result.getReference().length() > 1) {
            attrMap.put("END", result.getEnd());
            newRef = Allele.create(newRef.getBases()[0], true);
            gb.alleles(Arrays.asList(newRef, newRef));
        }
        //if GT is not homRef, correct it
        if (!isHomRefCall(result)) {
            gb.PL(new int[3]);  //3 for diploid PLs, automatically initializes to zero(?)
            gb.GQ(0);
            gb.alleles(Arrays.asList(newRef, newRef));
        }
        //there are some cases where there are low quality variants with homRef calls AND alt alleles!
        //TODO: the best thing would be to take the most likely alt's likelihoods
        else if (originalVC.getNAlleles() > 2) {
            final int[] idxVector = originalVC.getGLIndecesOfAlternateAllele(Allele.NON_REF_ALLELE);   //this is always length 3
            final int[] multiallelicPLs = g.getPL();
            final int[] newPLs = new int[3];
            newPLs[0] = multiallelicPLs[idxVector[0]];
            newPLs[1] = multiallelicPLs[idxVector[1]];
            newPLs[2] = multiallelicPLs[idxVector[2]];
            gb.PL(newPLs);
        }
        if (result.getGenotype(0).hasAD()) {
            int depth = (int) MathUtils.sum(result.getGenotype(0).getAD());
            gb.DP(depth);
            gb.attribute("MIN_DP", depth);
        }
        VariantContextBuilder builder = new VariantContextBuilder(result);
        return builder.alleles(Arrays.asList(newRef, Allele.NON_REF_ALLELE)).unfiltered().log10PError(VariantContext.NO_LOG10_PERROR).attributes(attrMap).genotypes(gb.make()).make(); //genotyping engine will add lowQual filter, so strip it off
    }

    private GenotypeBuilder makeGQ0RefCall(final VariantContext result, final Map<String, Object> attrMap) {
        Allele newRef = result.getReference();
        Genotype g = result.getGenotype(0);
        GenotypeBuilder gb = new GenotypeBuilder(g);
        //NB: If we're dropping a deletion allele, then we need to trim the reference and add an END tag with the vc stop position
        if (result.getReference().length() > 1) {
            attrMap.put("END", result.getEnd());
            newRef = Allele.create(newRef.getBases()[0], true);
            gb.alleles(Arrays.asList(newRef, newRef));
        }
        if (!isHomRefCall(result)) {
            gb.PL(new int[3]);  //3 for diploid PLs, automatically initializes to zero(?)
            gb.GQ(0).noAD().alleles(Arrays.asList(newRef, newRef)).noAttributes();
        }
        return gb;
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }
}
