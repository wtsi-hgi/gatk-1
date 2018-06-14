/*
* By downloading the PROGRAM you agree to the following terms of use:
*
* BROAD INSTITUTE
* SOFTWARE LICENSE AGREEMENT
* FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*
* This Agreement is made between the Broad Institute, Inc. with a principal address at 415 Main Street, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").
*
* WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
* WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
* NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*
* 1. DEFINITIONS
* 1.1 PROGRAM shall mean copyright in the object code and source code known as GATK3 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.org/gatk on the EFFECTIVE DATE.
*
* 2. LICENSE
* 2.1 Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
* The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only. For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
* 2.2 No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD. LICENSEE shall ensure that all of its users agree to the terms of this Agreement. LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
* 2.3 License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
*
* 3. PHONE-HOME FEATURE
* LICENSEE expressly acknowledges that the PROGRAM contains an embedded automatic reporting system ("PHONE-HOME") which is enabled by default upon download. Unless LICENSEE requests disablement of PHONE-HOME, LICENSEE agrees that BROAD may collect limited information transmitted by PHONE-HOME regarding LICENSEE and its use of the PROGRAM.  Such information shall include LICENSEE'S user identification, version number of the PROGRAM and tools being run, mode of analysis employed, and any error reports generated during run-time.  Collection of such information is used by BROAD solely to monitor usage rates, fulfill reporting requirements to BROAD funding agencies, drive improvements to the PROGRAM, and facilitate adjustments to PROGRAM-related documentation.
*
* 4. OWNERSHIP OF INTELLECTUAL PROPERTY
* LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies. LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
* Copyright 2012-2016 Broad Institute, Inc.
* Notice of attribution: The GATK3 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
* LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*
* 5. INDEMNIFICATION
* LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*
* 6. NO REPRESENTATIONS OR WARRANTIES
* THE PROGRAM IS DELIVERED AS IS. BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
* IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*
* 7. ASSIGNMENT
* This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*
* 8. MISCELLANEOUS
* 8.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
* 8.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
* 8.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
* 8.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested. All notices under this Agreement shall be deemed effective upon receipt.
* 8.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
* 8.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
* 8.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

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
 * java -jar GenomeAnalysisTK.jar \
 *   -T GenotypeGVCFs \
 *   -R reference.fasta \
 *   --variant sample1.g.vcf \
 *   --variant sample2.g.vcf \
 *   -o output.vcf
 * </pre>
 *
 * <h3>Caveats</h3>
 * <p>Only gVCF files produced by HaplotypeCaller (or CombineGVCFs) can be used as input for this tool. Some other
 * programs produce files that they call gVCFs but those lack some important information (accurate genotype likelihoods
 * for every position) that GenotypeGVCFs requires for its operation.</p>
 * <p>If the gVCF files contain allele specific annotations, add -G Standard -G AS_Standard to the command line.</p>
 *
 * <h3>Special note on ploidy</h3>
 * <p>This tool is able to handle any ploidy (or mix of ploidies) intelligently; there is no need to specify ploidy
 * for non-diploid organisms.</p>
 *
 */
public class ReblockGVCF extends VariantWalker {

    private static String GVCF_BLOCK = "GVCFBlock";
    private static int PLOIDY_TWO = 2;  //FIXME: hack to avoid extra arg collections

    /**
     * The gVCF file to reblock
     */
    @ArgumentCollection protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output(doc="File to which variants should be written")
    protected VariantContextWriter vcfWriter = null;

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

    /**
     * Which groups of annotations to add to the output VCF file. The single value 'none' removes the default group. See
     * the VariantAnnotator -list argument to view available groups. Note that this usage is not recommended because
     * it obscures the specific requirements of individual annotations. Any requirements that are not met (e.g. failing
     * to provide a pedigree file for a pedigree-based annotation) may cause the run to fail.
     */
    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls")
    protected List<String> annotationGroupsToUse = new ArrayList<>(Arrays.asList(new String[]{StandardAnnotation.class.getSimpleName()}));


    @Advanced
    @Argument(fullName="GVCFGQBands", shortName="GQB", doc="Exclusive upper bounds for reference confidence GQ bands " +
            "(must be in [1, 100] and specified in increasing order)")
    protected List<Integer> GVCFGQBands = new ArrayList<Integer>() {{
        add(20); add(100);
    }};

    @Advanced
    @Argument(fullName="drop_low_qual_variants", shortName="dropLowQuals", doc="Exclude variants that don't meet the reference GQ threshold from the reblocked GVCF to save space")
    protected boolean dropLowQuals = false;

    @Advanced
    @Argument(fullName="RGQ_threshold_to_drop", shortName="RGQthreshold", doc="Reference genotype quality (PL[0]) value below which to drop variant sites from the GVCF")
    protected double RGQthreshold = 0.0;

    @Advanced
    @Argument(fullName="do_qual_score_approximation", shortName="doQualApprox", doc="Add necessary INFO field annotation to perform QUAL approximation downstream")
    protected boolean doQualApprox = false;

    /**
     * The rsIDs from this file are used to populate the ID column of the output.  Also, the DB INFO flag will be set when appropriate. Note that dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    protected DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    // the genotyping engine
    private HaplotypeCallerGenotypingEngine genotypingEngine;
    // the INFO field annotation key names to remove
    private final List<String> infoFieldAnnotationKeyNamesToRemove = new ArrayList<>();

    // INFO Header names that require alt alleles
    final Set<String> infoHeaderAltAllelesLineNames = new LinkedHashSet<>();

    public void initialize() {
        final VCFHeader inputVCFHeader = getHeaderForVariants();
        final SampleList samples = new IndexedSampleList(inputVCFHeader.getGenotypeSamples()); //todo should this be getSampleNamesInOrder?
        final Set<VCFHeaderLine> headerLines = new LinkedHashSet<>(inputVCFHeader.getMetaDataInInputOrder());
        headerLines.addAll(getDefaultToolVCFHeaderLines());

        if ( dbsnp.dbsnp != null  ) {
            VCFStandardHeaderLines.addStandardInfoLines(headerLines, true, VCFConstants.DBSNP_KEY);
        }

        // create the genotyping engine
        //don't do allele-specific genotyping, i.e. don't output the AS_QUAL
        boolean doAlleleSpecificGenotyping = false;
        UnifiedArgumentCollection UAC = createUAC();
        UAC.outputMode = OutputMode.EMIT_ALL_SITES;
        UAC.annotateAllSitesWithPLs = true;
        genotypingEngine = new HaplotypeCallerGenotypingEngine(hcArgs, samplesList, FixedAFCalculatorProvider.createThreadSafeProvider(hcArgs), ! hcArgs.doNotRunPhysicalPhasing);
        genotypingEngine.setAnnotationEngine(annotationEngine);

        // take care of the VCF headers
        final Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), true);

        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));   // needed for gVCFs without DP tags
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.RAW_QUAL_APPROX_KEY));
        headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.VARIANT_DEPTH_KEY));
        headerLines.add(new VCFInfoHeaderLine("MQ_DP", 1, VCFHeaderLineType.Integer, "Depth over variant samples for better MQ calculation"));

        headerLines.remove(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.DOWNSAMPLED_KEY));
        headerLines.remove(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.HAPLOTYPE_SCORE_KEY));
        headerLines.remove(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.INBREEDING_COEFFICIENT_KEY));
        headerLines.remove(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_COUNT_KEY));
        headerLines.remove(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY));

        // Remove old GVCFBlock header entries
        for ( final Iterator<VCFHeaderLine> iter = headerLines.iterator(); iter.hasNext(); ) {
            if ( iter.next().getKey().contains(GVCF_BLOCK) ) {
                iter.remove();
            }
        }

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

        VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(
                new File(outputVCF),
                readsDictionary,
                createOutputVariantMD5,
                options.toArray(new Options[options.size()])
        );

        try {
            writer = new GVCFWriter(writer, hcArgs.GVCFGQBands, hcArgs.genotypeArgs.samplePloidy);
        } catch ( IllegalArgumentException e ) {
            throw new CommandLineException.BadArgumentValue("GQBands", "are malformed: " + e.getMessage());
        }        vcfWriter.writeHeader(vcfHeader);

        logger.info("Notice that the -ploidy parameter is ignored in " + getClass().getSimpleName() + " tool as this is automatically determined by the input variant files");
    }

    // get VariantContexts from input gVCFs, merge, and regenotype
    public void apply(VariantContext variant, ReadsContext reads, ReferenceContext ref, FeatureContext features) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return null;

        final List<VariantContext> vcsAtThisLocus = tracker.getValues(variantCollection.variants, context.getLocation());
        if ( vcsAtThisLocus == null || vcsAtThisLocus.isEmpty()) {
            return null;
        }

        final Byte refBase = INCLUDE_NON_VARIANTS ? ref.getBase() : null;
        final boolean removeNonRefSymbolicAllele = !INCLUDE_NON_VARIANTS;
        final VariantContext combinedVC = vcsAtThisLocus.get(0); //just take the first -- there should only be one if these are properly formed gVCFs (i.e. deletion positions aren't listed twice)
        return combinedVC == null ? null : regenotypeVC(tracker, ref, combinedVC);
    }

    /**
     * Re-genotype (and re-annotate) a combined genomic VC
     *
     * @param tracker        the ref tracker
     * @param ref            the ref context
     * @param originalVC     the combined genomic VC
     * @return a new VariantContext or null if the site turned monomorphic and we don't want such sites
     */
    protected VariantContext regenotypeVC(final RefMetaDataTracker tracker, final ReferenceContext ref, final VariantContext originalVC) {
        if ( originalVC == null ) {
            throw new IllegalArgumentException("originalVC cannot be null");
        } else if (!isProperlyPolymorphic(originalVC) && !INCLUDE_NON_VARIANTS) {
            return null;
        }

        VariantContext result = originalVC;

        //Pass back ref-conf homRef sites
        if (result.getLog10PError() == VariantContext.NO_LOG10_PERROR) {
            if (dropLowQuals && result.getGenotype(0).getGQ() == 0) {
                return null;
            }
            else if (result.getGenotype(0).isCalled() && result.getGenotype(0).isHomRef()) {

                return result;
            }
            else if (!result.getGenotype(0).isCalled() && result.getGenotype(0).hasPL() && result.getGenotype(0).getPL()[0] == 0) {
                return result;
            }
        }
        boolean isHomRefCall = result.getGenotype(0).isHomRef();

        //don't need to calculate quals for sites with no data whatsoever or sites already genotyped homRef
        if (result.getAttributeAsInt(VCFConstants.DEPTH_KEY,0) > 0 && !isHomRefCall) {
            result = genotypingEngine.calculateGenotypes(originalVC);
        }

        if (result == null || (!isProperlyPolymorphic(result) && !INCLUDE_NON_VARIANTS)) {
            return null;
        }

        //variants with PL[0] less than threshold get turned to homRef with PL=[0,0,0], shouldn't get INFO attributes
        //make sure we can call het variants with GQ >= RGQthreshold in joint calling downstream
        if(result.getGenotype(0).getPL()[0] < RGQthreshold) {
            if(dropLowQuals && !isHomRefCall) {
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
            if (!isHomRefCall) {
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
                int varDP = QualByDepth.getDepth(result.getGenotypes(), null, null);
                if (varDP == 0) {  //prevent QD=Infinity case
                    varDP = result.getAttributeAsInt(VCFConstants.DEPTH_KEY, Integer.MAX_VALUE);
                }
                attrMap.put(GATKVCFConstants.VARIANT_DEPTH_KEY, varDP);
            }
            VariantContextBuilder builder = new VariantContextBuilder(result);
            builder.attributes(attrMap);

            boolean allelesNeedSubsetting = false;
            List<Allele> allelesToDrop = new ArrayList<>();
            //only put in alleles that are called
            for(final Allele currAlt : result.getAlternateAlleles()) {
                boolean foundMatch = false;
                for(final Allele gtAllele: result.getGenotype(0).getAlleles()) {
                    if(gtAllele.equals(currAlt,false)) {
                        foundMatch = true;
                        break;
                    }
                    if(gtAllele.equals(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE)){
                        if(dropLowQuals) { //don't regenotype, just drop it -- this is a GQ 0 case if ever I saw one
                            return null;
                        }
                        else {
                            //TODO: extract this method to a "make GQ0 ref"
                            Allele newRef = result.getReference();
                            Genotype g = result.getGenotype(0);
                            GenotypeBuilder gb = new GenotypeBuilder(g);
                            //NB: If we're dropping a deletion allele, then we need to trim the reference and add an END tag with the vc stop position
                            if(result.getReference().length() > 1) {
                                attrMap.put("END", result.getEnd());
                                newRef = Allele.create(newRef.getBases()[0], true);
                                gb.alleles(Arrays.asList(newRef, newRef));
                            }
                            if (!isHomRefCall) {
                                gb.PL(new int[3]);  //3 for diploid PLs, automatically initializes to zero(?)
                                gb.GQ(0).noAD().alleles(Arrays.asList(newRef, newRef)).noAttributes();
                            }
                            return builder.alleles(Arrays.asList(newRef, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE)).unfiltered().log10PError(VariantContext.NO_LOG10_PERROR).attributes(attrMap).genotypes(gb.make()).make();
                        }
                    }
                }
                if(!foundMatch && !currAlt.isSymbolic()) {
                    allelesNeedSubsetting = true;
                    allelesToDrop.add(currAlt);
                }
            }
            //remove any AD reads for the non-ref
            int nonRefInd = result.getAlleleIndex(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE);
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
            attrMap.put("MQ_DP", originalVC.getAttributeAsInt(VCFConstants.DEPTH_KEY,0));

            if(allelesNeedSubsetting) {
                List<Allele> newAlleleSet = new ArrayList<>();
                for(final Allele a : result.getAlleles()) {
                    newAlleleSet.add(Allele.create(a,false));
                }
                newAlleleSet.removeAll(allelesToDrop);
                builder.alleles(newAlleleSet);
                if(!genotypesWereModified) {
                    builder.genotypes(GATKVariantContextUtils.subsetAlleles(result, newAlleleSet, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN));
                }
                else {  //again, only initialize a builder if we have to
                    builder.genotypes(GATKVariantContextUtils.subsetAlleles(new VariantContextBuilder(result).genotypes(newGenotypes).make(), newAlleleSet, GenotypeAssignmentMethod.USE_PLS_TO_ASSIGN));
                }
                return GATKVariantContextUtils.reverseTrimAlleles(builder.attributes(attrMap).make());
            }
            return builder.attributes(attrMap).genotypes(newGenotypes).make();
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

    /**
     * Creates a UnifiedArgumentCollection with appropriate values filled in from the arguments in this walker
     * @return a complete UnifiedArgumentCollection
     */
    private UnifiedArgumentCollection createUAC() {
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.genotypeArgs = genotypeArgs.clone();

        //whether to emit non-variant sites is not contained in genotypeArgs and must be passed to uac separately
        uac.outputMode = INCLUDE_NON_VARIANTS ? OutputMode.EMIT_ALL_CONFIDENT_SITES : OutputMode.EMIT_VARIANTS_ONLY;
        return uac;
    }

    public VariantContextWriter makeVCFWriter( final String outputVCF, final SAMSequenceDictionary readsDictionary,
                                               final boolean createOutputVariantIndex, final boolean  createOutputVariantMD5,
                                               final boolean sitesOnlyMode ) {
        Utils.nonNull(outputVCF);
        Utils.nonNull(readsDictionary);

        final List<Options> options = new ArrayList<>(2);
        if (createOutputVariantIndex) {options.add(Options.INDEX_ON_THE_FLY);}
        if (sitesOnlyMode) {options.add(Options.DO_NOT_WRITE_GENOTYPES);}

        VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(
                new File(outputVCF),
                readsDictionary,
                createOutputVariantMD5,
                options.toArray(new Options[options.size()])
        );

            try {
                writer = new GVCFWriter(writer, GVCFGQBands, genotypeArgs.samplePloidy);
            } catch ( IllegalArgumentException e ) {
                throw new CommandLineException.BadArgumentValue("GQBands", "are malformed: " + e.getMessage());
            }


        return writer;
    }
}
